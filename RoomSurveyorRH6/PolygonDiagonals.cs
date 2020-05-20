using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyorRH6
{
    public class PolygonDiagonals : GH_Component
    {
        public PolygonDiagonals()
          : base("Polygon Diagonals", "PolyDiag",
            "PolygonDiagonals description",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed planar polyline", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Internal Diagonals", "ID", "The internal diagonals of the polygon", GH_ParamAccess.list);
            pManager.AddLineParameter("External Diagonals", "ED", "The diagonals that are outside of the polygon", GH_ParamAccess.list);
            pManager.AddLineParameter("Intersecting Diagonals", "SID", "The diagonals of the polygon that intersect it", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();

            if (!DA.GetData(0, ref curve)) return;

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            if (curve.TryGetPolyline(out poly)) { }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P input expects a polyline to be provided");
                return;
            }

            if (poly.IsValid == false)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Closed polylines must have at least 2 segments");
                return;
            }

            if (poly.IsClosed == false)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }

            double[,] diagMatrix = RoomSurvey4.UpperMatrix(poly);


        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.PolyDiagonals_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("9c9b769b-a945-4b44-9955-e82cf5e75393"); }
        }
    }
}
