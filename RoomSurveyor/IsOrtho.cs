using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor

{
    public class IsOrtho : GH_Component
    {

        public IsOrtho()
          : base("Is OrthoPolygon", "IsOrtho",
            "Determines if a polygon or closed polyline is an orthopolygon, i.e. all its internal angles are either 90 or 270 degrees",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed polyline or a polygon", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBooleanParameter("Ortho", "O", "True if polygon is an orthopolygon, false if any of the polygon's internal angles is different from 90 or 270", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            bool ortho;

            if (!DA.GetData(0, ref curve)) return;

            if (curve.TryGetPolyline(out poly)) { }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P input expects a polygon to be provided");
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

            ortho = IsPolyOrtho(poly);

            DA.SetData(0, ortho);
        }

        /// <summary>
        /// Checks if all the angles in polyline are orthogonal, either 90 or 270.
        /// </summary>
        /// <returns><c>true</c>, if an ortho-polygon was input, <c>false</c> otherwise.</returns>
        /// <param name="polyline">Polyline.</param>
        public static bool IsPolyOrtho(Polyline polyline)
        {
            double turn = 0.0;
            bool ortho = false;

            for (int i = 0; i < polyline.SegmentCount; i++)
            {
                Point3d p0, p1, p2;
                if (i == 0)
                {
                    p0 = polyline[polyline.Count - 2];
                    p1 = polyline[i];
                    p2 = polyline[i + 1];
                }
                else
                {
                    p0 = polyline[i - 1];
                    p1 = polyline[i];
                    p2 = polyline[i + 1];
                }
                Vector3d v1 = p1 - p0;
                Vector3d v2 = p2 - p1;
                turn += Math.Abs(v1 * v2);
            }

            if (Math.Abs(turn) < 0.00001) { ortho = true; }

            return ortho;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.IsOrtho_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("b41bc25c-dc5b-44c3-8726-1555c3848a3f"); }
        }
    }
}
