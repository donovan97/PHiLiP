#include "delaunay.h"

namespace PHiLiP {
namespace ProperOrthogonalDecomposition {

Edge::Edge(RowVector2d A, RowVector2d B)
    : node1(A)
    , node2(B)
{}

bool Edge::operator==(const Edge &other) {
    return ((other.node1 == node1 && other.node2 == node2) ||
            (other.node1 == node2 && other.node2 == node1));
}

Triangle::Triangle(RowVector2d A, RowVector2d B, RowVector2d C)
    : node1(A)
    , node2(B)
    , node3(C)
    , edge1(A, B)
    , edge2(B, C)
    , edge3(A, C)
{
    //Formula for coordinates of the circumcenter
    double D = 2*(A(0)*(B(1) - C(1)) + B(0)*(C(1)-A(1)) + C(0)*(A(1)-B(1)));
    double center_x = (1/D)*((std::pow(A(0), 2) + std::pow(A(1), 2))*(B(1) - C(1)) + (std::pow(B(0), 2) + std::pow(B(1), 2))*(C(1) - A(1)) + (std::pow(C(0), 2) + std::pow(C(1), 2))*(A(1) - B(1)));
    double center_y = (1/D)*((std::pow(A(0), 2) + std::pow(A(1), 2))*(C(0) - B(0)) + (std::pow(B(0), 2) + std::pow(B(1), 2))*(A(0) - C(0)) + (std::pow(C(0), 2) + std::pow(C(1), 2))*(B(0) - A(0)));
    circumcenter(center_x, center_y);
    radius_squared = std::pow(A(0) - center_x, 2) + std::pow(A(1) - center_y, 2);
    edges = {edge1, edge2, edge3};
    nodes = {node1, node2, node3};
}

bool Triangle::operator==(const Triangle &other) {
    return ((other.node1 == node1 && other.node2 == node2 && other.node3 == node3) ||
            (other.node1 == node2 && other.node2 == node3 && other.node3 == node1) ||
            (other.node1 == node3 && other.node2 == node1 && other.node3 == node2));
}


Delaunay::Delaunay(Eigen::Matrix<double, Eigen::Dynamic, 2> points)
{
    triangulate(points);
}

void Delaunay::triangulate(Eigen::Matrix<double, Eigen::Dynamic, 2> points)
{
    RowVector2d min = points.colwise().minCoeff();
    RowVector2d max = points.colwise().maxCoeff();

    double dx = max(0) - min(0);
    double dy = max(0) - min(0);
    const auto dmax = std::max(dx, dy);
    const auto midx = (min(0) + max(0))/2;
    const auto midy = (min(1) + max(1))/2;

    RowVector2d super_node_1(midx - 20 * dmax, midy - dmax);
    RowVector2d super_node_2(midx, midy + 20 * dmax);
    RowVector2d super_node_3(midx + 20 * dmax, midy - dmax);

    triangulation.emplace_back(Triangle(super_node_1, super_node_2, super_node_3));

    for(auto point : points.rowwise()){
        std::vector<Triangle> bad_triangles;
        std::vector<Edge> bad_edges;

        for(const auto& tri : triangulation){
            //Check if point is inside the circumcircle
            double dist_squared = std::pow(tri.circumcenter(0) - point(0), 2) + std::pow(tri.circumcenter(1) - point(1), 2);
            if((dist_squared - tri.radius_squared) <= 1E-04){
                bad_triangles.push_back(tri);
                bad_edges.push_back(tri.edge1);
                bad_edges.push_back(tri.edge2);
                bad_edges.push_back(tri.edge3);
            }
        }

        std::vector<Edge> good_edges;

        for(auto& tri : bad_triangles){
            for(auto& tri_edge : tri.edges){
                for(auto& tri2: bad_triangles){
                    if ((tri == tri2)) {
                        continue;
                    }
                    for(auto& tri2_edge : tri2.edges){
                        if ((tri == tri2)) {
                            continue;
                        }
                    }
                    good_edges.push_back(tri_edge);
                }
            }
        }

        for(auto& bad_tri : bad_triangles) {
            //triangulation.erase(std::remove(triangulation.begin(), triangulation.end(), tri), triangulation.end());
            triangulation.erase(std::remove_if(triangulation.begin(), triangulation.end(),
                   [&bad_tri](Triangle& tri) {
                        return tri == bad_tri;
                    }), triangulation.end());
        }

        for(auto& edge : good_edges){
            Triangle tri(point, edge.node1, edge.node2);
            triangulation.push_back(tri);
        }
    }

    for(auto& tri : triangulation){
        for(auto& node : tri.nodes){
            if((node == super_node_1) || (node == super_node_2) || (node == super_node_3)){
                triangulation.erase(std::remove_if(triangulation.begin(), triangulation.end(),
                                                   [&tri](Triangle& tri_x) {
                                                       return tri == tri_x;
                                                   }), triangulation.end());
            }
        }
    }
}


}
}
