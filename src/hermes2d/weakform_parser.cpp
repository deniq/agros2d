// This file is part of Agros2D.
//
// Agros2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Agros2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Agros2D.  If not, see <http://www.gnu.org/licenses/>.
//
// hp-FEM group (http://hpfem.org/)
// University of Nevada, Reno (UNR) and University of West Bohemia, Pilsen
// Email: agros2d@googlegroups.com, home page: http://hpfem.org/agros2d/

#include "weakform_parser.h"
#include "scene.h"
#include "util.h"
#include "marker.h"
#include "hermes2d.h"
#include "module.h"
#include "coupling.h"
#include "module_agros.h"
#include "field.h"
#include "problem.h"
#include "datatable.h"

void ParserForm::initParserForm(FieldInfo *fieldInfo)
{
    m_fieldInfo = fieldInfo;
    m_couplingInfo = NULL;
    assert(fieldInfo);
    parser = new Parser(fieldInfo);
}

void ParserForm::initParserForm(CouplingInfo *couplingInfo)
{
    m_fieldInfo = NULL;
    m_couplingInfo = couplingInfo;
    assert(couplingInfo);
    parser = new Parser(couplingInfo);
}

ParserForm::~ParserForm()
{
    delete parser;
}


void ParserForm::initParser(QList<Material *> materials, Boundary *boundary)
{
    if(m_fieldInfo)
        parser->parser.push_back(m_fieldInfo->module()->expressionParser());
    else if(m_couplingInfo)
        parser->parser.push_back(m_couplingInfo->coupling()->expressionParser());
    else
        assert(0);

    // coordinates
    parser->parser[0]->DefineVar(Util::problem()->config()->labelX().toLower().toStdString(), &px);
    parser->parser[0]->DefineVar(Util::problem()->config()->labelY().toLower().toStdString(), &py);

    // current solution
    parser->parser[0]->DefineVar("uval", &puval);
    parser->parser[0]->DefineVar("ud" + Util::problem()->config()->labelX().toLower().toStdString() , &pudx);
    parser->parser[0]->DefineVar("ud" + Util::problem()->config()->labelY().toLower().toStdString(), &pudy);

    // test function
    parser->parser[0]->DefineVar("vval", &pvval);
    parser->parser[0]->DefineVar("vd" + Util::problem()->config()->labelX().toLower().toStdString(), &pvdx);
    parser->parser[0]->DefineVar("vd" + Util::problem()->config()->labelY().toLower().toStdString(), &pvdy);

    // previous solution
    parser->parser[0]->DefineVar("upval", &pupval);
    parser->parser[0]->DefineVar("upd" + Util::problem()->config()->labelX().toLower().toStdString(), &pupdx);
    parser->parser[0]->DefineVar("upd" + Util::problem()->config()->labelY().toLower().toStdString(), &pupdy);

    // solution from previous time level
    parser->parser[0]->DefineVar("uptval", &puptval);
    parser->parser[0]->DefineVar("uptd" + Util::problem()->config()->labelX().toLower().toStdString(), &puptdx);
    parser->parser[0]->DefineVar("uptd" + Util::problem()->config()->labelY().toLower().toStdString(), &puptdy);

    // time step
    parser->parser[0]->DefineVar("deltat", &pdeltat);

    // coupling
    for(int comp = 0; comp < maxSourceFieldComponents; comp++)
    {
        parser->parser[0]->DefineVar("source" + QString().setNum(comp).toStdString(), &source[comp]);
        parser->parser[0]->DefineVar("source" + QString().setNum(comp).toStdString() + "d" + Util::problem()->config()->labelX().toLower().toStdString(), &sourcedx[comp]);
        parser->parser[0]->DefineVar("source" + QString().setNum(comp).toStdString() + "d" + Util::problem()->config()->labelY().toLower().toStdString(), &sourcedy[comp]);
    }

    parser->setParserVariables(materials, boundary);

    for (std::map<QString, double>::iterator it = parser->parser_variables.begin(); it != parser->parser_variables.end(); ++it)
    {
        parser->parser[0]->DefineVar(it->first.toStdString(), &it->second);
    }
}

// **********************************************************************************************

template <typename Scalar>
CustomParserMatrixFormVol<Scalar>::CustomParserMatrixFormVol(unsigned int i, unsigned int j,
                                                             std::string area,
                                                             Hermes::Hermes2D::SymFlag sym,
                                                             std::string expression,
                                                             FieldInfo *fieldInfo,
                                                             CouplingInfo *couplingInfo,
                                                             Material *material1,
                                                             Material *material2)
//TODO kam vsude probubla material
// ->fieldInfo z materialu se v Parser form pouziva k projiti vsech fieldInfo->module->material_type_variables popr. module->get_boundary_type
// m_material .. pouzije se pro ziskani hodnot promennych
// initParser -> set parser variables ... take hodnoty promennych

    : Hermes::Hermes2D::MatrixFormVol<Scalar>(i, j, area, sym),
      m_fieldInfo(fieldInfo), m_couplingInfo(couplingInfo), m_material1(material1), m_material2(material2)
{
    if(fieldInfo)
        initParserForm(fieldInfo);
    else if(couplingInfo)
        initParserForm(couplingInfo);
    else
        assert(0);

    QList<Material *> materials;
    materials.push_back(material1);
    if(material2)
        materials.push_back(material2);
    initParser(materials, NULL);

    parser->parser[0]->SetExpr(expression);

    if(this->linearityType() == LinearityType_Linear)
    {
        pupval = 0;  // todo: ???
        if(m_fieldInfo)
        {
            foreach (Module::MaterialTypeVariable *variable, m_fieldInfo->module()->materialTypeVariables())
            {
                Value value = m_material1->getValue(variable->id);

                // table
                if (value.table->size() > 0)
                {
                    parser->parser_variables[variable->shortname] = value.value(pupval);
                    parser->parser_variables["d" + variable->shortname] = value.derivative(pupval);
                }

                // parser->parser_variables[variable->shortname] = m_material->get_value(variable->id).value(sqrt(pupdx*pupdx + pupdy*pupdy));
                // parser->parser_variables["d" + variable->shortname] = m_material->get_value(variable->id).derivative(sqrt(pupdx*pupdx + pupdy*pupdy));

                // if (variable->shortname == "mur")
                //     qDebug() << 1.0/parser->parser_variables[variable->shortname]/(4*M_PI*1e-7);
            }
        }
    }
}

template <typename Scalar>
Scalar CustomParserMatrixFormVol<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
                                                Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::ExtData<Scalar> *ext) const
{
    double result = 0;

    pdeltat = Util::problem()->config()->timeStep().number();

    for (int i = 0; i < n; i++)
    {
        px = e->x[i];
        py = e->y[i];

        puval = u->val[i];
        pudx = u->dx[i];
        pudy = u->dy[i];

        pvval = v->val[i];
        pvdx = v->dx[i];
        pvdy = v->dy[i];

        pupval = u_ext[this->i]->val[i];
        pupdx = u_ext[this->i]->dx[i];
        pupdy = u_ext[this->i]->dy[i];

        // previous solution
        if (this->linearityType() != LinearityType_Linear)
        {
            if(m_fieldInfo)
            {
                foreach (Module::MaterialTypeVariable *variable, m_fieldInfo->module()->materialTypeVariables())
                {
                    Value value = m_material1->getValue(variable->id);

                    // table
                    if (value.table->size() > 0)
                    {
                        parser->parser_variables[variable->shortname] = value.value(pupval);
                        parser->parser_variables["d" + variable->shortname] = value.derivative(pupval);
                    }

                    // parser->parser_variables[variable->shortname] = m_material->get_value(variable->id).value(sqrt(pupdx*pupdx + pupdy*pupdy));
                    // parser->parser_variables["d" + variable->shortname] = m_material->get_value(variable->id).derivative(sqrt(pupdx*pupdx + pupdy*pupdy));

                    // if (variable->shortname == "mur")
                    //     qDebug() << 1.0/parser->parser_variables[variable->shortname]/(4*M_PI*1e-7);
                }
            }
        }
        else
        {
            pupval = 0.0;
            pupdx = 0.0;
            pupdy = 0.0;
        }

        if(! m_material2)
        {
            if (m_fieldInfo && (m_fieldInfo->analysisType() == AnalysisType_Transient))
            {
                puptval = ext->fn[this->j]->val[i];
                puptdx = ext->fn[this->j]->dx[i];
                puptdy = ext->fn[this->j]->dy[i];
            }
        }

        try
        {
            result += wt[i] * parser->parser[0]->Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cout << "CustomParserMatrixFormVol: " << parser->parser[0]->GetExpr() << " - " << e.GetMsg() << std::endl;
        }
    }

    return result;
}

template <typename Scalar>
Hermes::Ord CustomParserMatrixFormVol<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
                                                   Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::ExtData<Hermes::Ord> *ext) const
{
    return Hermes::Ord(20);
    Hermes::Ord result;

    for (int i = 0; i < n; i++)
        result += u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i];

    return result;
}

template <typename Scalar>
CustomParserMatrixFormVol<Scalar>* CustomParserMatrixFormVol<Scalar>::clone()
{
    return new CustomParserMatrixFormVol(this->i, this->j, this->areas[0], (Hermes::Hermes2D::SymFlag) this->sym, parser->parser[0]->GetExpr(),
                                         this->m_fieldInfo, this->m_couplingInfo, this->m_material1, this->m_material2);
}

template <typename Scalar>
LinearityType CustomParserMatrixFormVol<Scalar>::linearityType() const
{
    LinearityType linearityType;

    if(m_fieldInfo)
        linearityType = m_fieldInfo->linearityType();
    else
        linearityType = m_couplingInfo->linearityType();

    return linearityType;
}


template <typename Scalar>
CustomParserVectorFormVol<Scalar>::CustomParserVectorFormVol(unsigned int i, unsigned int j,
                                                             std::string area,
                                                             std::string expression,
                                                             FieldInfo *fieldInfo,
                                                             CouplingInfo *couplingInfo,
                                                             Material *material1,
                                                             Material *material2)
    : Hermes::Hermes2D::VectorFormVol<Scalar>(i, area), m_fieldInfo(fieldInfo), m_couplingInfo(couplingInfo),
      m_material1(material1), m_material2(material2), j(j)
{
    if(fieldInfo)
        initParserForm(fieldInfo);
    else if(couplingInfo)
        initParserForm(couplingInfo);
    else
        assert(0);

    QList<Material *> materials;
    materials.push_back(material1);
    if (material2)
        materials.append(material2);
    initParser(materials, NULL);

    parser->parser[0]->SetExpr(expression);

    if(this->linearityType() == LinearityType_Linear)
    {
        pupval = 0;  // todo: ???
        if(m_fieldInfo)
        {
            foreach (Module::MaterialTypeVariable *variable, m_fieldInfo->module()->materialTypeVariables())
            {
                Value value = m_material1->getValue(variable->id);

                // table
                if (value.table->size() > 0)
                {
                    parser->parser_variables[variable->shortname] = m_material1->getValue(variable->id).value(pupval);
                }

                // parser->parser_variables[variable->shortname] = m_material->get_value(variable->id).value(sqrt(pupdx*pupdx + pupdy*pupdy));

                // if (variable->shortname == "epsr")
                //     qDebug() << parser->parser_variables[variable->shortname];
            }
        }
    }
}

template <typename Scalar>
Scalar CustomParserVectorFormVol<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
                                                Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::ExtData<Scalar> *ext) const
{
    double result = 0;

    pdeltat = Util::problem()->config()->timeStep().number();

    for (int i = 0; i < n; i++)
    {
        px = e->x[i];
        py = e->y[i];

        pvval = v->val[i];
        pvdx = v->dx[i];
        pvdy = v->dy[i];

        pupval = u_ext[this->j]->val[i];
        pupdx = u_ext[this->j]->dx[i];
        pupdy = u_ext[this->j]->dy[i];

        if (this->linearityType() != LinearityType_Linear)
        {
            if(m_fieldInfo)
            {
                foreach (Module::MaterialTypeVariable *variable, m_fieldInfo->module()->materialTypeVariables())
                {
                    Value value = m_material1->getValue(variable->id);

                    // table
                    if (value.table->size() > 0)
                    {
                        parser->parser_variables[variable->shortname] = m_material1->getValue(variable->id).value(pupval);
                    }

                    // parser->parser_variables[variable->shortname] = m_material->get_value(variable->id).value(sqrt(pupdx*pupdx + pupdy*pupdy));

                    // if (variable->shortname == "epsr")
                    //     qDebug() << parser->parser_variables[variable->shortname];
                }
            }
        }

        if(m_material2){
            // we have material2 -> it is coupling form, without time parameter
            if (ext->get_nf() > 0)
            {
                for(int comp = 0; comp < ext->get_nf(); comp++)
                {
                    source[comp] = ext->fn[comp]->val[i];
                    sourcedx[comp] = ext->fn[comp]->dx[i];
                    sourcedy[comp] = ext->fn[comp]->dy[i];
                }
            }
        }
        else
        {
            if (m_fieldInfo && (m_fieldInfo->analysisType() == AnalysisType_Transient))
            {
                puptval = ext->fn[this->j]->val[i];
                puptdx = ext->fn[this->j]->dx[i];
                puptdy = ext->fn[this->j]->dy[i];
            }
        }

        try
        {
            result += wt[i] * parser->parser[0]->Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cout << "CustomParserVectorFormVol: " << parser->parser[0]->GetExpr() << " - " << e.GetMsg() << std::endl;
        }
    }

    return result;
}

template <typename Scalar>
Hermes::Ord CustomParserVectorFormVol<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
                                                   Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::ExtData<Hermes::Ord> *ext) const
{
    return Hermes::Ord(20);
    Hermes::Ord result;

    for (int i = 0; i < n; i++)
        result += v->dx[i] + v->dy[i];

    return result;
}


template <typename Scalar>
CustomParserVectorFormVol<Scalar>* CustomParserVectorFormVol<Scalar>::clone()
{
    return new CustomParserVectorFormVol(this->i, this->j, this->areas[0], parser->parser[0]->GetExpr(),
                                         this->m_fieldInfo, this->m_couplingInfo, this->m_material1, this->m_material2);
}

template <typename Scalar>
LinearityType CustomParserVectorFormVol<Scalar>::linearityType() const
{
    LinearityType linearityType;

    if(m_fieldInfo)
        linearityType = m_fieldInfo->linearityType();
    else
        linearityType = m_couplingInfo->linearityType();

    return linearityType;
}



// **********************************************************************************************

template <typename Scalar>
CustomParserMatrixFormSurf<Scalar>::CustomParserMatrixFormSurf(unsigned int i, unsigned int j,
                                                               std::string area, std::string expression,
                                                               Boundary *boundary)
    : Hermes::Hermes2D::MatrixFormSurf<Scalar>(i, j, area), m_boundary(boundary)
{
    initParserForm(boundary->getFieldInfo());
    initParser(QList<Material *>(), boundary);

    parser->parser[0]->SetExpr(expression);
}

template <typename Scalar>
Scalar CustomParserMatrixFormSurf<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v,
                                                 Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::ExtData<Scalar> *ext) const
{
    double result = 0;

    pdeltat = Util::problem()->config()->timeStep().number();

    for (int i = 0; i < n; i++)
    {
        px = e->x[i];
        py = e->y[i];

        puval = u->val[i];
        pudx = u->dx[i];
        pudy = u->dy[i];

        pvval = v->val[i];
        pvdx = v->dx[i];
        pvdy = v->dy[i];

        // todo: now we use only Newton solver
        // previous solution
//        if (m_fieldInfo->linearityType != LinearityType_Linear)
//        {
            pupval = u_ext[this->j]->val[i];
            pupdx = u_ext[this->j]->dx[i];
            pupdy = u_ext[this->j]->dy[i];
//        }
//        else
//        {
//            pupval = 0.0;
//            pupdx = 0.0;
//            pupdy = 0.0;
//        }

        try
        {
            result += wt[i] * parser->parser[0]->Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cout << "CustomParserMatrixFormSurf: " << parser->parser[0]->GetExpr() << " - " << e.GetMsg() << std::endl;
        }
    }

    return result;
}

template <typename Scalar>
Hermes::Ord CustomParserMatrixFormSurf<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v,
                                                    Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::ExtData<Hermes::Ord> *ext) const
{
    return Hermes::Ord(20);
    Hermes::Ord result;

    for (int i = 0; i < n; i++)
        result += u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i];

    return result;
}

template <typename Scalar>
CustomParserMatrixFormSurf<Scalar>* CustomParserMatrixFormSurf<Scalar>::clone()
{
    return new CustomParserMatrixFormSurf(this->i, this->j, this->areas[0], parser->parser[0]->GetExpr(),
                                          this->m_boundary);
}

template <typename Scalar>
CustomParserVectorFormSurf<Scalar>::CustomParserVectorFormSurf(unsigned int i, unsigned int j,
                                                               std::string area, std::string expression,
                                                               Boundary *boundary)
    : Hermes::Hermes2D::VectorFormSurf<Scalar>(i, area), j(j), m_boundary(boundary)
{
    initParserForm(boundary->getFieldInfo());
    initParser(QList<Material *>(), boundary);

    parser->parser[0]->SetExpr(expression);
}

template <typename Scalar>
Scalar CustomParserVectorFormSurf<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
                                                 Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::ExtData<Scalar> *ext) const
{
    double result = 0;

    pdeltat = Util::problem()->config()->timeStep().number();

    for (int i = 0; i < n; i++)
    {
        px = e->x[i];
        py = e->y[i];

        pvval = v->val[i];
        pvdx = v->dx[i];
        pvdy = v->dy[i];

        // todo: now we use only Newton solver
//        // previous solution
//        if (m_fieldInfo->linearityType != LinearityType_Linear)
//        {
            pupval = u_ext[this->j]->val[i];
            pupdx = u_ext[this->j]->dx[i];
            pupdy = u_ext[this->j]->dy[i];
//        }
//        else
//        {
//            pupval = 0.0;
//            pupdx = 0.0;
//            pupdy = 0.0;
//        }

        try
        {
//            if(parser->parser[0]->Eval())
//            {
//            cout << "area " <<  this->areas[0] << "  ";
//                for(std::map<std::string, double>::iterator it = parser->parser_variables.begin(); it != parser->parser_variables.end(); ++it)
//                {
//                    cout << (*it).first << "=" << (*it).second << ", ";
//                }
//                cout << "*****expresion***** " << parser->parser[0]->GetExpr() << ", value -> " << parser->parser[0]->Eval() << endl;
//            }

            result += wt[i] * parser->parser[0]->Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cout << "CustomParserVectorFormSurf: " << parser->parser[0]->GetExpr() << " - " << e.GetMsg() << std::endl;
        }
    }

    return result;
}

template <typename Scalar>
Hermes::Ord CustomParserVectorFormSurf<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
                                                    Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::ExtData<Hermes::Ord> *ext) const
{
    return Hermes::Ord(20);
    Hermes::Ord result;

    for (int i = 0; i < n; i++)
        result += v->dx[i] + v->dy[i];

    return result;
}

template <typename Scalar>
CustomParserVectorFormSurf<Scalar>* CustomParserVectorFormSurf<Scalar>::clone()
{
    return new CustomParserVectorFormSurf(this->i, this->j, this->areas[0], parser->parser[0]->GetExpr(),
                                          this->m_boundary);
}

// **********************************************************************************************

template <typename Scalar>
CustomExactSolution<Scalar>::CustomExactSolution(Hermes::Hermes2D::Mesh *mesh, std::string expression, Boundary *boundary)
    : Hermes::Hermes2D::ExactSolutionScalar<Scalar>(mesh)
{
    initParserForm(boundary->getFieldInfo());
    initParser(QList<Material *>(), boundary);

    parser->parser[0]->SetExpr(expression);
}

template <typename Scalar>
Scalar CustomExactSolution<Scalar>::value(double x, double y) const
{
    double result = 0;

    px = x;
    py = y;

    try
    {
        result = parser->parser[0]->Eval();
    }
    catch (mu::Parser::exception_type &e)
    {
        std::cout << "CustomExactSolution: " << parser->parser[0]->GetExpr() << " - " << e.GetMsg() << std::endl;
    }

    return result;
}

template <typename Scalar>
void CustomExactSolution<Scalar>::derivatives (double x, double y, Scalar& dx, Scalar& dy) const
{
}

template class CustomParserMatrixFormVol<double>;
template class CustomParserMatrixFormSurf<double>;
template class CustomParserVectorFormVol<double>;
template class CustomParserVectorFormSurf<double>;
template class CustomExactSolution<double>;
