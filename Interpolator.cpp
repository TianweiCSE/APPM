#include "Interpolator.h"

extern std::string input_dir;

Interpolator::Interpolator() {

}

Interpolator::Interpolator(const PrimalMesh* primal, const DualMesh* dual)
: primal(primal), dual(dual) {
    initInterpolateTensor_E();
    initInterpolateTensor_B();
    // test();
}

Interpolator::~Interpolator() {

}

const Tensor3& Interpolator::get_E_interpolator() const {
    return E_interpolator;
}

const Tensor3& Interpolator::get_B_interpolator() const {
    return B_interpolator;
}

void Interpolator::initInterpolateTensor_E() {

    // !!Note: Different from definition in the report, the size of the first dim of interpolation tensor is extended
    //         from the number of dual fluid face to the number of all dual faces. This is for the convenience of 
    //         implementation, since the discrete faraday laws (4.32) (4.33) are defined on each dual face.
    E_interpolator = Tensor3(dual->getCells().size(), dual->getNumberOfFaces());
    for (const Cell* cell : dual->getCells()) {
        const int nFaces = cell->getFaceList().size();
        Eigen::MatrixX3d mat(nFaces, 3);
        mat.setZero();
        for (int i = 0; i < nFaces; i++) {
            const Face* face = cell->getFaceList()[i];
            if (face->isBoundary()) {  // The associated variable is d^\partial, namely, dual face flux of D-field.
                if (face->isPlane()) { // plane face
                    mat.row(i) = face->getNormal() * face->getArea();
                }
                else {  // non-plane face
                    for (const Face* subf : face->getSubFaceList()) {
                        mat.row(i) += subf->getNormal() * subf->getArea();
                    }
                }
            }
            else { // The associated variable is e, namely, primal edge integral of E-field. 
                    // Note: dual face idx = primal edge idx
                const Edge* edge = primal->getEdge(face->getIndex());
                mat.row(i) = edge->getDirection() * edge->getLength();
            }
        }
        Eigen::Matrix3Xd temp = (mat.transpose() * mat).inverse() * mat.transpose();
        for (int i = 0; i < nFaces; i++) {
            const Face* face = cell->getFaceList()[i];
            E_interpolator.insert(cell->getIndex(), face->getIndex(), temp.col(i));
        }
    }
    E_interpolator = E_interpolator;
}

void Interpolator::initInterpolateTensor_B() {

    // !!Note: Different from definition in the report, the size of the first dim of interpolation tensor is extended
    //         from the number of dual fluid face to the number of all dual faces. This is for the convenience of 
    //         implementation, since the discrete faraday laws (4.32) (4.33) are defined on each dual face.  
    B_interpolator = Tensor3(dual->getCells().size(), dual->getNumberOfEdges());
    for (const Cell* cell : dual->getCells()) {
        const int nEdges = cell->getEdgeList().size();
        Eigen::MatrixX3d mat(nEdges, 3);
        for (int i = 0; i < nEdges; i++) {
            const Edge* edge = cell->getEdgeList()[i];
            if (edge->isBoundary()) { // The associated variable is h^\partial, namely, dual edge intergral of H-field.
                mat.row(i) = edge->getDirection() * edge->getProjectedLength();
            }
            else { // The associated variable is b, namely, primal face flux of B-field.
                    // Note: dual edge index = primal face index
                const Face* face = primal->getFace(edge->getIndex());
                assert(face->isPlane());
                mat.row(i) = face->getNormal() * face->getArea();
            }
        }
        Eigen::Matrix3Xd temp = (mat.transpose() * mat).inverse() * mat.transpose();
        for (int i = 0; i < nEdges; i++) {
            const Edge* edge = cell->getEdgeList()[i];
            B_interpolator.insert(cell->getIndex(), edge->getIndex(), temp.col(i));
        }
    }
}

void Interpolator::test() const {
    { // Test E_interpolator
        Eigen::Vector3d E_vec;
        E_vec << 1.0, 2.0, 3.0;
        Eigen::VectorXd e(primal->getNumberOfEdges());
        Eigen::VectorXd dp(dual->facet_counts.nF_boundary);
        dp.setZero();
        for (int i = 0; i < dual->getNumberOfFaces(); i++) {
            const Face* face = dual->getFace(i);
            if (face->isBoundary()) {
                if (face->isPlane()) {
                    dp[i - primal->getNumberOfEdges()] = face->getNormal().dot(E_vec) * face->getArea();
                }
                else {
                    for (const Face* subf : face->getSubFaceList()) {
                        dp[i - primal->getNumberOfEdges()] += subf->getNormal().dot(E_vec) * subf->getArea();
                    }
                }
            }
            else {
                const Edge* edge = primal->getEdge(i);
                e[i] = edge->getDirection().dot(E_vec) * edge->getLength();
            }
        }
        Eigen::VectorXd temp(e.size() + dp.size());
        temp << e, dp; // concatenate [e, dp]^T
        Eigen::MatrixXd E_int = get_E_interpolator().oneContract(temp);
        std::ofstream file( "test_E_int.dat");
        file << E_int;
    }
    { // Test B_interpolator
        Eigen::Vector3d B_vec;
        B_vec << 1.0, 2.0, 3.0;
        Eigen::VectorXd b(primal->getNumberOfFaces());
        Eigen::VectorXd hp(dual->facet_counts.nE_boundary);
        for (int i = 0; i < dual->getNumberOfEdges(); i++) {
            const Edge* edge = dual->getEdge(i);
            if (edge->isBoundary()) {
                hp[i - primal->getNumberOfFaces()] = edge->getDirection().dot(B_vec) * edge->getProjectedLength();
            }
            else {
                const Face* face = primal->getFace(i);
                b[i] = face->getNormal().dot(B_vec) * face->getArea();
            }
        }
        Eigen::VectorXd temp(b.size() + hp.size());
	    temp << b, hp;  // concatenate [b, hp]^T
	    Eigen::MatrixXd B_int = get_B_interpolator().oneContract(temp);
        std::ofstream file( "test_B_int.dat");
        file << B_int;
    }

}