#pragma once

#include "GeometryItem.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"

class Edge2: public Edge
{
public:
    Edge2();
    Edge2(const Edge* e1, const Edge* e2);
    ~Edge2();

    // Return the directions of both edges 
    const Eigen::MatrixXd getDirection() const;
    // Return the sum of two lengths
    const double getLength() const;

    bool hasConnectedFace(const std::vector<Edge*> & faceEdges) const; 
    // Return the joint point
    const Eigen::Vector3d getHalfwayPosition() const;


private:
    const Edge* _e1;
    const Edge* _e2;

};




