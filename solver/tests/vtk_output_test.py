"""Compare AS3's two visualization formats using the VTK readers."""

import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def run_solver(executable, example, folder, method, output_format, fields, precision, start_time, time_step):
    """Write initial, scheduled, and final snapshots in a temporary directory."""
    text = (example / "param.cfg").read_text()
    grid = re.search(r"^GRID_DIRECTORY\s*=\s*(\S+)", text, re.MULTILINE)[1].rstrip("/")
    options = {
        "TYPE_DG_METHOD": method,
        "TYPE_QUADRATURE": "LEGENDRE_GAUSS_LOBATTO" if method == "DG_SEM" else "GAUSS_LEGENDRE",
        "BONUS_QUADRATURE_FACTOR": 0,
        "OUTPUT_VIS_FORMAT": output_format,
        "OUTPUT_VIS_FILENAME": "vis/flow",
        "OUTPUT_SOL_FILENAME": "sol/flow",
        "WRITE_VIS_VARIABLE": fields,
        "START_TIME": start_time,
        "TIME_STEP": time_step,
        "WRITE_VIS_FREQ": 2,
        "MAX_ITER_TIME": 3,
    }
    for name, value in options.items():
        text = re.sub(rf"^[ \t]*{name}[ \t]*=.*\n?", "", text, flags=re.MULTILINE)
        text += f"\n{name} = {value}\n"
    text = re.sub(r"^[ \t]*OUTPUT_VIS_PRECISION[ \t]*=.*\n?", "", text, flags=re.MULTILINE)
    if precision is not None:
        text += f"\nOUTPUT_VIS_PRECISION = {precision}\n"

    (folder / "vis").mkdir(parents=True)
    (folder / "sol").mkdir()
    (folder / grid).symlink_to(example / grid, target_is_directory=True)
    (folder / "param.cfg").write_text(text)
    result = subprocess.run(
        [str(executable), "param.cfg"], cwd=folder, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        env=dict(os.environ, OMP_NUM_THREADS="2", OMP_DYNAMIC="FALSE"), timeout=90,
    )
    if result.returncode:
        raise RuntimeError(result.stdout)


def read_grid(path, binary_xml, expected_time):
    """Read the grid and check its physical time."""
    reader = vtkXMLUnstructuredGridReader() if binary_xml else vtkUnstructuredGridReader()
    if not binary_xml:
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
    reader.SetFileName(str(path))
    if binary_xml:
        reader.UpdateInformation()
        info = reader.GetOutputInformation(0)
        key = vtkStreamingDemandDrivenPipeline.TIME_STEPS()
        assert info.Length(key) == 1 and info.Get(key, 0) == expected_time, path
    reader.Update()
    grid = reader.GetOutput()
    assert reader.GetErrorCode() == 0 and grid.GetNumberOfCells() > 0, path
    time_array = grid.GetFieldData().GetArray("TimeValue")
    assert time_array is not None, path
    time_values = vtk_to_numpy(time_array)
    assert time_values.dtype == np.float64, path
    np.testing.assert_array_equal(time_values, [expected_time], err_msg=str(path))
    return grid


def compare_grids(legacy, xml, expected, dtype):
    """Check coordinates, cell topology, and every requested point-field tuple."""
    old_points, new_points = (vtk_to_numpy(grid.GetPoints().GetData()) for grid in (legacy, xml))
    assert old_points.dtype == new_points.dtype == dtype
    assert new_points.shape == (xml.GetNumberOfPoints(), 3) and not np.any(new_points[:, 2])
    np.testing.assert_array_equal(new_points, old_points)
    for getter in ("GetConnectivityArray", "GetOffsetsArray"):
        np.testing.assert_array_equal(vtk_to_numpy(getattr(xml.GetCells(), getter)()),
                                      vtk_to_numpy(getattr(legacy.GetCells(), getter)()))
    assert xml.GetNumberOfCells() == legacy.GetNumberOfCells()
    assert all(xml.GetCellType(i) == legacy.GetCellType(i) for i in range(xml.GetNumberOfCells()))

    a, b = legacy.GetPointData(), xml.GetPointData()
    for fields in (a, b):
        names = {fields.GetArrayName(i) for i in range(fields.GetNumberOfArrays())}
        assert names == expected, (names, expected)
    for name in expected:
        old, new = (vtk_to_numpy(fields.GetArray(name)) for fields in (a, b))
        assert old.dtype == new.dtype == dtype and np.isfinite(new).all(), name
        np.testing.assert_array_equal(new, old, err_msg=name)


def check_series_indexes(root, snapshot_times):
    """Check that the indexes map physical times to relative snapshot paths."""
    legacy = root / "VTK_LEGACY_BINARY" / "vis" / "flow.vtk.series"
    xml = root / "VTK_XML_BINARY" / "vis" / "flow.pvd"
    series = json.loads(legacy.read_text())
    assert series["file-series-version"] == "1.0"
    assert series["files"] == [
        {"name": f"flow_{index}.vtk", "time": time} for index, time in enumerate(snapshot_times)
    ]
    datasets = ET.parse(xml).findall("./Collection/DataSet")
    assert [(data.attrib["file"], float(data.attrib["timestep"])) for data in datasets] == [
        (f"flow_{index}.vtu", time) for index, time in enumerate(snapshot_times)
    ]


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: vtk_output_test.py /path/to/AS3")
    executable = Path(sys.argv[1]).resolve()
    examples = Path(__file__).resolve().parents[2] / "example"
    cases = (
        ("NODAL_DG", "( DENSITY, MOMENTUM, TOTAL_ENERGY, PRESSURE, VELOCITY, VORTICITY, MACH, TEMPERATURE, ENTROPY )",
         {"Density", "Momentum", "Energy", "Pressure", "Velocity", "Vorticity", "Mach", "Temperature", "Entropy"}),
        ("DG_SEM", "( VELOCITY, ENTROPY, MOMENTUM, PRESSURE )", {"Velocity", "Entropy", "Momentum", "Pressure"}),
    )
    start_time = 0.125
    time_step = 5.0e-5
    second_step = start_time + time_step + time_step
    snapshot_times = (start_time, second_step, second_step + time_step)
    with tempfile.TemporaryDirectory(prefix="as3-vtk-test-") as temp:
        for method, fields, expected in cases:
            default_grids = []
            for precision, dtype in ((None, np.float32), ("SINGLE", np.float32), ("DOUBLE", np.float64)):
                root = Path(temp) / method / (precision or "default")
                for format_name in ("VTK_LEGACY_BINARY", "VTK_XML_BINARY"):
                    run_solver(executable, examples / "curved_p5", root / format_name,
                               method, format_name, fields, precision, start_time, time_step)
                check_series_indexes(root, snapshot_times)
                for index, expected_time in enumerate(snapshot_times):
                    legacy = read_grid(root / "VTK_LEGACY_BINARY" / "vis" / f"flow_{index}.vtk",
                                       False, expected_time)
                    xml = read_grid(root / "VTK_XML_BINARY" / "vis" / f"flow_{index}.vtu",
                                    True, expected_time)
                    compare_grids(legacy, xml, expected, dtype)
                    if precision is None:
                        default_grids.append(xml)
                    else:
                        reference = default_grids[index]
                        points = vtk_to_numpy(xml.GetPoints().GetData())
                        rounded = points.astype(np.float32)
                        np.testing.assert_array_equal(rounded, vtk_to_numpy(reference.GetPoints().GetData()))
                        if precision == "DOUBLE":
                            assert np.any(points != rounded), "Coordinates were rounded before double output."
                        retains_precision = False
                        for name in expected:
                            data = vtk_to_numpy(xml.GetPointData().GetArray(name))
                            old = vtk_to_numpy(reference.GetPointData().GetArray(name))
                            rounded = data.astype(np.float32)
                            np.testing.assert_array_equal(rounded, old, err_msg=name)
                            retains_precision |= bool(np.any(data != rounded))
                        if precision == "DOUBLE":
                            assert retains_precision, "Point fields were rounded before double output."
                print(f"{method}, {precision or 'default SINGLE'}: all snapshots and times match across formats.")


if __name__ == "__main__":
    main()
