#include "geometry_structure.hpp"

#include <array>
#include <iostream>
#include <stdexcept>

namespace
{
void Require(bool condition, const char *message)
{
  if (!condition)
    throw std::runtime_error(message);
}

bool SameAddress(CFlattenedFaceIndex a, CFlattenedFaceIndex b)
{
  return a.mFaceType == b.mFaceType && a.mIndexBatch == b.mIndexBatch &&
         a.mIndexFace == b.mIndexFace;
}

void CheckEmptyBatches()
{
  struct Batch
  {
    std::size_t count;
    std::size_t GetnFace() const { return count; }
  };
  CFaceCollection<Batch> faces;
  Require(faces.GetnFace() == 0 && faces.GetnBatch() == 0, "Default collection is not empty");
  faces.Initialize({{0}, {2}, {0}, {3}, {0}});
  Require(faces.GetnFace() == 5 && faces.GetnBatch() == 5, "Empty batches changed counts");
  Require(faces.ResolveFace(0).mIndexBatch == 1, "Leading empty batch was selected");
  Require(faces.ResolveFace(1).mIndexFace == 1, "First batch endpoint is wrong");
  Require(faces.ResolveFace(2).mIndexBatch == 3, "Interior empty batch was selected");
  Require(faces.ResolveFace(4).mIndexFace == 2, "Last batch endpoint is wrong");
  faces.Initialize({{1}});
  Require(faces.GetnFace() == 1 && faces.GetnBatch() == 1, "Old ranges survived reinitialization");
  faces.Initialize({});
  Require(faces.GetnFace() == 0 && faces.GetnBatch() == 0, "Reinitialization did not clear ranges");
}

void CheckThinZones()
{
  for (const auto dims : {std::array<std::size_t, 2>{1, 1}, {1, 4}, {4, 1}})
  {
    CSinglezoneGeometry zone("test", 0, 1, dims[0], dims[1], ETypeDOF::LGL, true);
    CInteriorFaceBatch iFaces(ETypeFace::IFACE, &zone);
    CInteriorFaceBatch jFaces(ETypeFace::JFACE, &zone);
    Require(iFaces.GetnFace() == dims[1] * (dims[0] - 1), "Thin-zone i-face count is wrong");
    Require(jFaces.GetnFace() == dims[0] * (dims[1] - 1), "Thin-zone j-face count is wrong");
  }
}

void CheckTopology(const CMultizoneGeometry &geometry)
{
  // Every local edge of every element must occur exactly once in the topology.
  // Two-sided faces consume two local edges; physical boundaries consume one.
  as3vector2d<std::array<unsigned int, 4>> uses(geometry.GetnZone());
  for (std::size_t zone = 0; zone < uses.size(); ++zone)
    uses[zone].resize(geometry.GetSinglezoneGeometry(zone)->GetnElem(), {0, 0, 0, 0});

  const auto mark = [&](std::size_t zone, std::size_t element, EFaceLocation location)
  {
    Require(zone < uses.size(), "Face references a missing zone");
    Require(element < uses[zone].size(), "Face references a missing element");
    Require(++uses[zone][element][static_cast<std::size_t>(location)] == 1,
            "An element edge occurs in more than one face");
  };

  for (int direction = 0; direction < 2; ++direction)
  {
    const auto &faces =
        direction == 0 ? geometry.GetMultizoneIFaces() : geometry.GetMultizoneJFaces();
    const auto maxFace = direction == 0 ? EFaceLocation::IMAX : EFaceLocation::JMAX;
    const auto minFace = direction == 0 ? EFaceLocation::IMIN : EFaceLocation::JMIN;
    for (const auto &batch : faces.GetInteriorFaces().GetBatches())
      for (const auto &face : batch.GetFaces())
      {
        mark(batch.GetZoneIndex(), face.mIndexElementMinus, maxFace);
        mark(batch.GetZoneIndex(), face.mIndexElementMinus + batch.GetNeighborElementOffset(),
             minFace);
      }
    for (const auto &batch : faces.GetBoundaryFaces().GetBatches())
      for (const auto &face : batch.GetFaces())
        mark(batch.GetZoneIndex(), face.mIndexElement, batch.GetFaceLocation());
    for (const auto &batch : faces.GetInterfaceFaces().GetBatches())
      for (const auto &face : batch.GetFaces())
      {
        Require(batch.GetOwnerFaceDirection() == faces.GetFaceDirection(),
                "Interface is stored in the wrong logical direction");
        mark(batch.GetOwnerZoneIndex(), face.mIndexElementOwner, batch.GetOwnerFaceLocation());
        mark(batch.GetNeighborZoneIndex(), face.mIndexElementNeighbor,
             batch.GetNeighborFaceLocation());
      }

    // Check the flat order against independent enumeration of the batches.
    std::size_t index = 0;
    const auto checkCollection = [&](const auto &collection, ETypeFaceGeometry kind)
    {
      for (std::size_t batch = 0; batch < collection.GetnBatch(); ++batch)
        for (std::size_t face = 0; face < collection.GetBatch(batch).GetnFace(); ++face)
        {
          const CFlattenedFaceIndex expected{batch, face, kind};
          const auto cached = direction == 0 ? geometry.GetFlattenedIndexIFace(index)
                                             : geometry.GetFlattenedIndexJFace(index);
          Require(SameAddress(cached, expected), "Flattened face order differs from batch order");
          ++index;
        }
    };
    checkCollection(faces.GetInteriorFaces(), ETypeFaceGeometry::INTERNAL);
    checkCollection(faces.GetBoundaryFaces(), ETypeFaceGeometry::BOUNDARY);
    checkCollection(faces.GetInterfaceFaces(), ETypeFaceGeometry::INTERFACE);
    Require(index == faces.GetnFacesTotal(), "Flattened face count is wrong");

    const auto &permutation = direction == 0 ? geometry.GetIFaceLoadBalancedPermutation()
                                             : geometry.GetJFaceLoadBalancedPermutation();
    as3vector1d<bool> seen(index, false);
    for (std::size_t i = 0; i < index; ++i)
    {
      const auto mapped = permutation.GetPermutationIndex(i);
      Require(mapped < index && !seen[mapped], "Load balancing lost or duplicated a face");
      seen[mapped] = true;
    }

    int mismatches = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static) reduction(+ : mismatches)
#endif
    for (std::size_t i = 0; i < index; ++i)
    {
      const auto cached = direction == 0 ? geometry.GetFlattenedIndexIFaceLoadBalanced(i)
                                         : geometry.GetFlattenedIndexJFaceLoadBalanced(i);
      if (!SameAddress(cached, faces.ResolveFace(permutation.GetPermutationIndex(i))))
        ++mismatches;
    }
    Require(mismatches == 0, "Parallel face lookup differs from cached addresses");
  }
  for (const auto &zone : uses)
    for (const auto &element : zone)
      for (auto count : element)
        Require(count == 1, "An element edge is missing from the topology");
}
} // namespace

int main(int argc, char **argv)
{
  if (argc != 2)
    return 2;
  try
  {
    CheckEmptyBatches();
    CheckThinZones();
    CConfig config(argv[1]);
    COpenMP openmp(&config);
    CMultizoneGeometry geometry(&config);
    geometry.InitializeGridTopology(&openmp);
    CheckTopology(geometry);
    std::cout << "Face connectivity checks passed with " << openmp.GetnThread() << " thread(s).\n";
  }
  catch (const std::exception &error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
