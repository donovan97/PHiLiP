#ifndef __DELAUNAY__
#define __DELAUNAY__

#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include <algorithm>

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVector2d;

/// Delaunay triangulation
class Edge
{
public:
    /// Constructor
    Edge(RowVector2d node1, RowVector2d node2);

    /// Destructor
    ~Edge() {};

    RowVector2d node1;
    RowVector2d node2;
    bool operator==(const Edge& other);
};

class Triangle
{
public:
    /// Constructor
    Triangle(RowVector2d node1, RowVector2d node2, RowVector2d node3);

    /// Destructor
    ~Triangle() {};

    RowVector2d circumcenter;
    RowVector2d node1;
    RowVector2d node2;
    RowVector2d node3;
    Edge edge1;
    Edge edge2;
    Edge edge3;
    double radius_squared;
    std::vector<Edge> edges;
    std::vector<RowVector2d> nodes;

    bool operator==(const Triangle& other);
};

/// Delaunay triangulation
class Delaunay
{
public:
    /// Constructor
    Delaunay(Eigen::Matrix<double, Eigen::Dynamic, 2> points);

    /// Destructor
    ~Delaunay() {};

    void triangulate(Eigen::Matrix<double, Eigen::Dynamic, 2> points);

    std::vector<Triangle> triangulation;

};

}
}

#endif
