#include "Interpolator.h"

Interpolator::Interpolator() {

}

Interpolator::Interpolator(const PrimalMesh* primal, const DualMesh* dual)
: primal(primal), dual(dual) {
    initInterpolateTensor_E();
    initInterpolateTensor_B();
}

Interpolator::~Interpolator() {

}

const Tensor3& Interpolator::get_E_interpolator() const {
    return *E_interpolator;
}

const Tensor3& Interpolator::get_B_interpolator() const {
    return *B_interpolator;
}

void Interpolator::initInterpolateTensor_E() {

    // !!Note: Different from definition in the report, the size of the first dim of interpolation tensor is extended
    //         from the number of dual fluid face to the number of all dual faces. This is for the convenience of 
    //         implementation, since the discrete faraday laws (4.32) (4.33) are defined on each dual face.
    E_interpolator = new Tensor3(dual->getCells().size(), dual->getNumberOfFaces());
    for (const Cell* cell : dual->getFluidCells()) {
        double area_sum = 0;
        for (const Face* face : cell->getFaceList()) {
            // In the first traverse, collect sum of face areas. This is for computing the weight of each face
            area_sum += face->getArea();
        }
        for (const Face* face : cell->getFaceList()) {
            // In the second traverse, construct the tuple3 entries
            if (face->isBoundary()) {  // The associated variable is d^\partial, namely, dual face flux of D-field.
                E_interpolator->insert(cell->getIndex(), face->getIndex(), face->getNormal() / area_sum);
            }
            else { // The associated variable is e, namely, primal edge integral of E-field. 
                    // Note: dual face idx = primal edge idx
                E_interpolator->insert(
                    cell->getIndex(), 
                    face->getIndex(), 
                    face->getNormal() / primal->getEdge(face->getIndex())->getLength() * face->getArea() / area_sum
                );
            }
        }
    }
}

void Interpolator::initInterpolateTensor_B() {

    // !!Note: Different from definition in the report, the size of the first dim of interpolation tensor is extended
    //         from the number of dual fluid face to the number of all dual faces. This is for the convenience of 
    //         implementation, since the discrete faraday laws (4.32) (4.33) are defined on each dual face.  
    B_interpolator = new Tensor3(dual->getCells().size(), dual->getNumberOfEdges());
    for (const Cell* cell : dual->getFluidCells()) {
        double length_sum = 0;
        for (const Edge* edge : cell->getEdgeList()) {
            // In the first travese, collect the sum of edge lengths
            length_sum += edge->getLength();
        }
        for (const Edge* edge : cell->getEdgeList()) {
            // In the second traverse, construct the tuple3 entries
            if (edge->isBoundary()) { // The associated variable is h^\partial, namely, dual edge intergral of H-field.
                B_interpolator->insert(cell->getIndex(), edge->getIndex(), edge->getDirection() / length_sum);
            }
            else { // The associated variable is b, namely, primal face flux of B-field.
                    // Note: dual edge index = primal face index
                B_interpolator->insert(
                    cell->getIndex(),
                    edge->getIndex(),
                    edge->getDirection() / primal->getFace(edge->getIndex())->getArea() * edge->getLength() / length_sum
                );
            }
        }
    }
}