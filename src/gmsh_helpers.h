#pragma once

#include <filesystem>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "helpers.h"

#include "gmsh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

void parseMsh(
  const std::filesystem::path mshPath, 
  std::unique_ptr<ManifoldSurfaceMesh>& mesh, 
  std::unique_ptr<VertexPositionGeometry>& geometry, 
  globalBoundaryConditions& globalBdyConditions);