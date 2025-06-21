#include "Smoothing.hh"
#include "./eigen3/Eigen/Dense"
#include <vector>
using VertexVertexIter = TriMesh::VertexVertexIter;
using HalfedgeHandle   = TriMesh::HalfedgeHandle;
using EdgeHandle       = TriMesh::EdgeHandle;




// ======================================================================
// EXERCISE 1.1
// ========================================================================
void Smoothing::uniform_smooth(const int _iterations) {
    //Point laplacian;
    //int w = 0.;
    std::vector<Point> points_new(mesh_.n_vertices(), Point(0.));
    for (auto vh : mesh_.vertices())
        points_new[vh.idx()] = mesh_.point(vh);

    for (int iter=0; iter<_iterations; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        //w=1 so I do not change the initialized w and neither use it, just ignoring it because the weights are 1.
        for (auto vh : mesh_.vertices()){
             Point laplacian(0.0, 0.0, 0.0);
             if (!mesh_.is_boundary(vh)){
                 int neighbor_count=0;
                 for (VertexVertexIter vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
                    laplacian += mesh_.point(*vv_it) - mesh_.point(vh);
                    neighbor_count += 1;
                }
                 
                 if (neighbor_count > 0) {
                      laplacian /= neighbor_count;
                      // assuming lambda = 1
                      points_new[vh.idx()] = mesh_.point(vh) + laplacian;
                  } else {
                      points_new[vh.idx()] = mesh_.point(vh);  // keep unchanged
                  }
              }
            }
        }
       
        // ------------- IMPLEMENT HERE ---------
        // update vertex positions
        
        for (auto vh : mesh_.vertices()) {
            mesh_.point(vh) = points_new[vh.idx()];
        }
    

    mesh_.update_normals();
}

// ======================================================================
// EXERCISE 1.2
// ========================================================================
void Smoothing::cotan_laplacian_smooth(const int _iterations) {
    Point laplacian;
    double ww = 0.;
    std::vector<Point> points_new(mesh_.n_vertices(), Point(0.));
    for (auto vh : mesh_.vertices())
        points_new[vh.idx()] = mesh_.point(vh);

    for (int iter=0; iter<_iterations; ++iter) {
        // update edge weights
        calc_edges_weights();
        // ------------- IMPLEMENT HERE ---------
        for (auto vh : mesh_.vertices()){
        double sum_wij{0};
        Point laplacian(0.0, 0.0, 0.0);
            if (!mesh_.is_boundary(vh)){
                for (VertexVertexIter vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
                    HalfedgeHandle heh = mesh_.find_halfedge(vh, *vv_it);
                    EdgeHandle eh = mesh_.edge_handle(heh);
                    double wij = mesh_.property(edge_weight_, eh); 

                    laplacian += wij*(mesh_.point(*vv_it) - mesh_.point(vh));
                    sum_wij+=wij; 
                }
                if(sum_wij>1e-8) {  // safety check
                    laplacian /= sum_wij;
                    points_new[vh.idx()] = mesh_.point(vh) + laplacian;
                }
            }
        }   
        
    }
        

        // ------------- IMPLEMENT HERE ---------
        // update vertex positions
        for (auto vh : mesh_.vertices()) {
            mesh_.point(vh) = points_new[vh.idx()];
        }
    

    mesh_.update_normals();
}

// ======================================================================
// EXERCISE 2
// ========================================================================
void Smoothing::implicit_smooth(const double _timestep) {
    int n = (int)mesh_.n_vertices();
    std::vector<Point> points_new(mesh_.n_vertices(), Point(0.0));
    // compute cotan edge weights and vertex areas
    calc_weights();
    // ------------- IMPLEMENT HERE ---------
    //store Ai
    std::vector<double> A(n, 0.0); 

    for (auto vh : mesh_.vertices()) {
         double area = mesh_.property(vertex_weight_, vh);
         A[vh.idx()] = area ;
    }
    //build D inverse
    Eigen::VectorXd D(n);
    for (int i = 0; i < n; ++i)
        D(i) =  1/(2.0 * A[i]);
    // Prepare to construct M
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    for (auto vh : mesh_.vertices()){
        int i=vh.idx();
        double diag=0.0;
        if (!mesh_.is_boundary(vh)){
            for (VertexVertexIter vv_it = mesh_.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
                  HalfedgeHandle heh = mesh_.find_halfedge(vh, *vv_it);
                  EdgeHandle eh = mesh_.edge_handle(heh);
                  int j = vv_it->idx();
                  double wij = mesh_.property(edge_weight_, eh); 
                  diag-=wij;
                  triplets.push_back(T(i,j,wij));
            }
            triplets.push_back(T(i,i,diag));
        }
    }
    Eigen::SparseMatrix<double> M(n, n);
    M.setFromTriplets(triplets.begin(), triplets.end());
  
    for (auto vh : mesh_.vertices())
        points_new[vh.idx()] = mesh_.point(vh);
    //(D.inverse()-timestep*lambda_ M)points_new=D.inverse()*points_previous
    //convert points_previous to a matrix
    Eigen::MatrixXd points_previous(mesh_.n_vertices(),3);
    for(auto vh: mesh_.vertices()){
        int i=vh.idx();
        Point p = mesh_.point(vh);
        points_previous.row(i) = Eigen::Vector3d(p[0], p[1], p[2]);
    }
    // Convert D (a dense vector) into a sparse diagonal matrix
    Eigen::SparseMatrix<double> Dinv_sparse(n, n);
    std::vector<Eigen::Triplet<double>> Dinv_triplets;
    for (int i = 0; i < n; ++i)
        Dinv_triplets.emplace_back(i, i, 1.0 / D(i));
    Dinv_sparse.setFromTriplets(Dinv_triplets.begin(), Dinv_triplets.end());

    //build right hand side term D.inverse()*points_previous
    Eigen::MatrixXd RHS= Dinv_sparse*points_previous;

    //build left hand side (D.inverse()-timestep*M)points_new
    Eigen::SparseMatrix<double> LHS = Dinv_sparse - _timestep * M;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(LHS);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver failed!" << std::endl;
        return;
    }

    //  Solve for each dimension separately
    Eigen::MatrixXd Points_new(mesh_.n_vertices(), 3);
    for (int dim = 0; dim < 3; ++dim)
        Points_new.col(dim) = solver.solve(RHS.col(dim));

    //  Copy back to mesh
    for (auto vh : mesh_.vertices()) {
        int i = vh.idx();
        mesh_.set_point(vh, Point(Points_new(i, 0), Points_new(i, 1), Points_new(i, 2)));
    }
    
    // ------------- IMPLEMENT HERE ---------
    mesh_.update_normals();
}

// ======================================================================
// EXERCISE 3.1
// ========================================================================
void Smoothing::uniform_laplacian_enhance_feature(const int _iterations, const int _coefficient) {
    // ------------- IMPLEMENT HERE ---------
    std::vector<Point> points_new(mesh_.n_vertices(), Point(0.));
    std::vector<Point> points_old(mesh_.n_vertices(), Point(0.));
    std::vector<Point> difference(mesh_.n_vertices(), Point(0.));
    std::vector<Point> points_enhanced(mesh_.n_vertices(), Point(0.));
    //store the old points in vector points_old
    for (auto vh : mesh_.vertices())
        points_old[vh.idx()]=mesh_.point(vh);
    //smooth the points
    uniform_smooth(_iterations);
    //store the new points in vector points_new
    for (auto vh : mesh_.vertices()){
        points_new[vh.idx()]=mesh_.point(vh);
        difference[vh.idx()]=points_new[vh.idx()]-points_old[vh.idx()];
        points_enhanced[vh.idx()]=points_old[vh.idx()]-_coefficient*difference[vh.idx()];
    }
    
    

    // ------------- IMPLEMENT HERE ---------
    for (auto vh : mesh_.vertices()) {
            mesh_.point(vh) = points_enhanced[vh.idx()];
        }
    mesh_.update_normals();
    
}

// ======================================================================
// EXERCISE 3.2
// ========================================================================
void Smoothing::cotan_laplacian_enhance_feature(const int _iterations, const int _coefficient) {
    // ------------- IMPLEMENT HERE ---------
    std::vector<Point> points_new(mesh_.n_vertices(), Point(0.));
    std::vector<Point> points_old(mesh_.n_vertices(), Point(0.));
    std::vector<Point> difference(mesh_.n_vertices(), Point(0.));
    std::vector<Point> points_enhanced(mesh_.n_vertices(), Point(0.));
    //store the old points in vector points_old
    for (auto vh : mesh_.vertices())
        points_old[vh.idx()]=mesh_.point(vh);
    //smooth the points
    cotan_laplacian_smooth(_iterations);
    //store the new points in vector points_new
    for (auto vh : mesh_.vertices()){
        points_new[vh.idx()]=mesh_.point(vh);
        difference[vh.idx()]=points_new[vh.idx()]-points_old[vh.idx()];
        points_enhanced[vh.idx()]=points_old[vh.idx()]-_coefficient*difference[vh.idx()];
    }

    // ------------- IMPLEMENT HERE ---------
    for (auto vh : mesh_.vertices()) {
            mesh_.point(vh) = points_enhanced[vh.idx()];
        }
    mesh_.update_normals();
}
