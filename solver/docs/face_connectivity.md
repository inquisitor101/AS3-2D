# Face connectivity {#face_connectivity}

The face data describe connectivity. Coordinates, normals, and metric terms live
with the element geometry. `face_structure.hpp` keeps the existing project naming
style while giving each level a specific meaning.

| Type | Responsibility |
| --- | --- |
| `CInteriorFace`, `CBoundaryFace`, `CInterfaceFace` | Element indices for one face. |
| `CInteriorFaceBatch` | All interior faces in one zone and one logical direction. |
| `CBoundaryFaceBatch` | All faces on one physical boundary marker: one zone edge and one boundary name. |
| `CInterfaceFaceBatch` | Paired faces on one connected pair of zone edges, including periodic connections. |
| `CFaceCollection<TBatch>` | Owns batches of one kind and maps a flat face index to its batch and local face. |
| `CDirectionalFaceConnectivity` | Owns interior, boundary, and interface collections across all zones for one logical direction. |

A batch shares metadata across its faces. For example, zone IDs and edge locations
are stored once for an interface batch, while each interface face stores only the
two element indices. A batch is independent of OpenMP thread assignments and
load-balancing chunks.

The two interface sides are **owner** and **neighbor**. These correspond to the
input marker's I and J sides. Logical **i/j** always refers to a grid direction.
The owner determines which directional container holds an interface; the neighbor
can have a different direction in a rotated zone. Owner/neighbor also works for
periodic connections between two edges of the same zone.

Within a zone, an interior face retains the existing minus/plus convention. The
face stores the minus element index. Its plus element is obtained by adding the
batch's neighbor element offset: 1 for i-faces, `niElem` for j-faces.

There are three distinct index spaces:

1. A direction-local index visits `[interior][boundary][interface]` in that order.
2. A collection-local index visits all faces of one kind across its batches.
3. A batch-local index selects one entry in that batch's face vector.

`CFaceCollection::ResolveFace` maps (2) to `(batch, local face)` in `CFaceAddress`.
`CDirectionalFaceConnectivity::ResolveFace` maps (1) to
`(kind, batch, local face)` in `CFlattenedFaceIndex`. The existing
`ETypeFaceGeometry::INTERNAL` value denotes interior faces.

Category offsets are private to the directional resolver. Geometry caches the
resolved addresses during topology initialization, before the parallel residual
loops. The existing load-balancing permutation indexes that cache. No batch search
was added to the residual loops.

Typical interface access is now:

```cpp
const auto address = geometry.GetFlattenedIndexIFaceLoadBalanced(i);
const auto& interfaces = geometry.GetMultizoneIFaces().GetInterfaceFaces();
const auto& batch = interfaces.GetBatch(address.mIndexBatch);
const auto& face = batch.GetFace(address.mIndexFace);

const auto ownerZone = batch.GetOwnerZoneIndex();
const auto ownerElement = face.mIndexElementOwner;
```

This example assumes `address.mFaceType == ETypeFaceGeometry::INTERFACE`, as in
the existing residual dispatch. A batch index is local to its face kind and its
directional container.

Collections expose const batch access. Changing their contents requires
`Initialize`, which rebuilds the flat index ranges. Initialize topology before
parallel traversal; the existing solver residual storage and synchronization
rules still apply.

The main API renames are:

| Previous name | Current name |
| --- | --- |
| `CMultizoneFaceGeometry` | `CDirectionalFaceConnectivity` |
| `CInternalFacesFamily` | `CInteriorFaceBatch` |
| `CBoundaryFacesFamily` | `CBoundaryFaceBatch` |
| `CInterfaceFacesFamily` | `CInterfaceFaceBatch` |
| `CGroupFaces<TFamily>` | `CFaceCollection<TBatch>` |
| `CElementFaceIndex` | `CFaceAddress` |
| `mIndexFamily` in face addresses | `mIndexBatch` |
| `Get*FacesFamily(index)` | `Get*Faces().GetBatch(index)` |
| `GetiZone()` / `GetjZone()` on interfaces | `GetOwnerZoneIndex()` / `GetNeighborZoneIndex()` |

Input marker families retain their names: they group input markers and can
produce multiple face batches. Their format and configuration syntax are
unchanged. `InitializeGridTopology` now takes only the OpenMP container; face
construction gets its marker information from the geometry.

Run `bash tests/run_face_connectivity.sh` from the project root to check the
connectivity independently of the full solver. The default input is
`example/curved_p5/param.cfg`. Use
`bash tests/run_face_connectivity.sh example/orthogonal_square_p5` for the
orthogonal grid.
The test runner uses Homebrew Clang and libomp by default; `CXX` and `OMP_PREFIX`
can select other installations on macOS. It builds a temporary executable with
debug checks, then tests the curved nine-zone grid with periodic, mixed, and
physical boundaries, each with one and four OpenMP threads. It also checks empty
batches, collection reinitialization, and zones one element wide.

The release OpenMP build succeeds with `make -f Makefile_Clang -j4`. The updated
boundary kernels, residual signatures, and accumulation paths are preserved.
`CBoundaryFaceBatch::IsMaxFace()` exposes the cached boundary flag used to choose
and accumulate the residual buffer.

End-to-end regression checks passed against a fresh build of the pre-refactor
sources. They cover the curved periodic example for 200 steps,
the curved wall and inlet/outlet examples for five steps, and the orthogonal
periodic and inlet/outlet examples for five steps, with one and four threads.

All 26 VTK frames matched byte-for-byte, and every numerical monitor row matched
for each configuration and thread count.
