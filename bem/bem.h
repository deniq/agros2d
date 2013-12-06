#ifndef BEM_H
#define BEM_H

#include <QList>
#include <QString>
#include <QObject>

#include "../agros2d-library/meshgenerator.h"
#include "../agros2d-library/hermes2d/field.h"
#include "../agros2d-library/hermes2d/problem.h"
#include "../agros2d-library/hermes2d/problem_config.h"
#include "../hermes2d/include/function/exact_solution.h"
#include "../agros2d-library/hermes2d/solver.h"
#include "bem_matrix.h"
#include "mesh.h"


class Bem
{

public:
    Bem(FieldInfo* field, std::tr1::shared_ptr<Hermes::Hermes2D::Mesh> mesh);
    ~Bem() {}
    void readMesh();
    void addPhysics();
    void assemblyMatrices();
    void solve();
    void domainSolution();
    double getValue(double x, double y);

    double quad(int oder, int j, Node node, Segment segment, double (Bem::*kernel)(Node, Segment, double));
    double gaussLaguerre(int oder, QList<Node> nodes, Segment segment,double (Bem::*kernel)(Node, Segment, double));

    double kernel_length(Node refNode, Segment segment, double xi);
    double kernel_laplace2D_derivation(Node refNode, Segment segment, double xi);
    double kernel_laplace2D(Node refNode, Segment segment, double xi);

    QString toString();
    Node integral(Node v, Node a, Node b);
    double potential(double x, double y);
    Node globalCoordinates(double xi, Segment segment);
    BemVector shapeFunction(int n, double xi);
    BemVector shapeFunction2D(int polyOrder, double s, double t);
    BemVector shapeFunctionDerivative(int polyOrder, double xi);
    Node normalVector(double xi, Segment segment);
    double jacobian(int polyOrder, double xi, Segment segment);

    // Private ?
    Mesh mesh;
    unsigned int m_polyOrder;

private:
    MeshSharedPtr m_mesh;
    FieldInfo * m_fieldInfo;
};

template <typename Scalar>
class BemSolution: public Hermes::Hermes2D::ExactSolutionScalar<Scalar>
{
public:
    BemSolution<Scalar>(MeshSharedPtr mesh);
    virtual ~BemSolution() {}
    virtual Scalar value (double x, double y) const;
    virtual void derivatives (double x, double y, Scalar& dx, Scalar& dy) const;
    virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const { return Hermes::Ord(20); }
    virtual Hermes::Hermes2D::MeshFunction<Scalar>* clone() const;    
    /// Saves the exact solution to an XML file.
    void save(const char* filename) const { }    
    void setSolver(Bem * bem) { m_bem = bem; }    

protected:    
    // virtual void precalculate(int order, int mask) { qDebug() << "OK";}
    double constant;
    Bem * m_bem;
};




#endif // BEM_H
