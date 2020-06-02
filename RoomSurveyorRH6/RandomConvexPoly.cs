using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyorr

{
    public class RandomConvexPoly : GH_Component
    {
        public RandomConvexPoly()
          : base("RandomConvexPoly", "RCP",
                 "Generates a random convex polygon. A C# version of the algorithm implemented in Java by Sander Verdonschot, based on Pavel Valtr's proof of the probability of convexity of a random point set.",
            "RoomSurveyor", "Polygon Generators")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Number", "N", "An integer with the number of corners of the polygon to generate", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "The random polygon", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int corners = 0;
            List<Point2d> pts = new List<Point2d>();
            List<Point3d> pts3D = new List<Point3d>();
            Polyline poly = new Polyline();

            if (!DA.GetData(0, ref corners)) return;

            if (corners < 3)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "N must be larger or equal to 3");
                return;
            }
            pts.AddRange(ValtrAlgorithm.GenerateRandomConvexPolygon(corners));

            foreach (Point2d p in pts)
            {
                Point3d pt3D = new Point3d(p.X, p.Y, 0);
                pts3D.Add(pt3D);
            }

            poly.AddRange(pts3D);
            poly.Add(pts3D[0]);

            DA.SetData(0, poly);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.RandomConvexPoly_Icon;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("36dfcb8b-a8b8-43cc-835b-44156a11bb2c"); }
        }
    }
}
