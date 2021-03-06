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

#include "{{ID}}_filter.h"
#include "{{ID}}_interface.h"

#include "util.h"
#include "util/global.h"

#include "hermes2d/problem.h"
#include "hermes2d/problem_config.h"

#include "hermes2d/plugin_interface.h"

{{CLASS}}ViewScalarFilter::{{CLASS}}ViewScalarFilter(const FieldInfo *fieldInfo, int timeStep, int adaptivityStep, SolutionMode solutionType,
                                           Hermes::vector<Hermes::Hermes2D::MeshFunctionSharedPtr<double> > sln,
                                           const QString &variable,
                                           PhysicFieldVariableComp physicFieldVariableComp)
    : Hermes::Hermes2D::Filter<double>(sln), m_fieldInfo(fieldInfo), m_timeStep(timeStep), m_adaptivityStep(adaptivityStep), m_solutionType(solutionType),
      m_variable(variable), m_physicFieldVariableComp(physicFieldVariableComp)
{
    m_variableHash = qHash(m_variable);

    {{#SPECIAL_FUNCTION_SOURCE}}
    if(m_fieldInfo->functionUsedInAnalysis("{{SPECIAL_FUNCTION_ID}}"))
        {{SPECIAL_FUNCTION_NAME}} = QSharedPointer<{{SPECIAL_EXT_FUNCTION_FULL_NAME}}>(new {{SPECIAL_EXT_FUNCTION_FULL_NAME}}(m_fieldInfo, 0));
    {{/SPECIAL_FUNCTION_SOURCE}}

    value = new double*[this->num];
    dudx = new double*[this->num];
    dudy = new double*[this->num];

    m_coordinateType = Agros2D::problem()->config()->coordinateType();
    m_labels = Agros2D::scene()->labels;
}

{{CLASS}}ViewScalarFilter::~{{CLASS}}ViewScalarFilter()
{
    delete [] value;
    delete [] dudx;
    delete [] dudy;
}

Hermes::Hermes2D::Func<double> *{{CLASS}}ViewScalarFilter::get_pt_value(double x, double y, bool use_MeshHashGrid, Hermes::Hermes2D::Element* e)
{
    return NULL;
}

void {{CLASS}}ViewScalarFilter::precalculate(int order, int mask)
{
    Hermes::Hermes2D::Quad2D* quad = this->quads[Hermes::Hermes2D::Function<double>::cur_quad];
    int np = quad->get_num_points(order, this->get_active_element()->get_mode());
    
    for (int k = 0; k < this->num; k++)
    {
        this->sln[k]->set_quad_order(order, Hermes::Hermes2D::H2D_FN_DEFAULT);
        /// \todo Find out why Qt & OpenGL renders the outputs color-less if dudx, dudy, valus are 'const double*'.
        dudx[k] = const_cast<double*>(this->sln[k]->get_dx_values());
        dudy[k] = const_cast<double*>(this->sln[k]->get_dy_values());
        value[k] = const_cast<double*>(this->sln[k]->get_fn_values());
    }

    this->update_refmap();

    double *x = this->refmap->get_phys_x(order);
    double *y = this->refmap->get_phys_y(order);
    Hermes::Hermes2D::Element *e = this->refmap->get_active_element();

    // set material
    SceneMaterial *material = m_labels->at(atoi(m_fieldInfo->initialMesh()->get_element_markers_conversion().
                                           get_user_marker(e->marker).marker.c_str()))->marker(m_fieldInfo);

    int elementMarker = e->marker;

    {{#VARIABLE_MATERIAL}}const Value *material_{{MATERIAL_VARIABLE}} = material->valueNakedPtr(QLatin1String("{{MATERIAL_VARIABLE}}"));
    {{/VARIABLE_MATERIAL}}    
    {{#VARIABLE_SOURCE}}
    if ((m_variableHash == {{VARIABLE_HASH}})
            && (m_coordinateType == {{COORDINATE_TYPE}})
            && (m_fieldInfo->analysisType() == {{ANALYSIS_TYPE}})
            && (m_physicFieldVariableComp == {{PHYSICFIELDVARIABLECOMP_TYPE}}))
        for (int i = 0; i < np; i++)
            this->values[0][0][i] = {{EXPRESSION}};
    {{/VARIABLE_SOURCE}}
}

{{CLASS}}ViewScalarFilter* {{CLASS}}ViewScalarFilter::clone() const
{
    Hermes::vector<Hermes::Hermes2D::MeshFunctionSharedPtr<double> > slns;

    for (int i = 0; i < this->num; i++)
        slns.push_back(this->sln[i]->clone());

    {{CLASS}}ViewScalarFilter *filter = new {{CLASS}}ViewScalarFilter(m_fieldInfo, m_timeStep, m_adaptivityStep, m_solutionType, slns, m_variable, m_physicFieldVariableComp);

    return filter;
}

