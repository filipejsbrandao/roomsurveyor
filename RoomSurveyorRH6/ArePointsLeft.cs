using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyorRH6
{
    public class ArePointsLeft : GH_Component
    {

        public ArePointsLeft()
          : base("ArePointsLeft", "AreLeft",
            "Determines if the points are left, on or right of the infinite line from P0 to P1 in 2D space",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("First Point", "P0", "The origin point", GH_ParamAccess.item);
            pManager.AddPointParameter("Second Point", "P1", "The second point determines the direction of the line", GH_ParamAccess.item);
            pManager.AddPointParameter("Points", "Pt", "Points to test", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Left", "L", "Points that are left of the line from P0 to P1", GH_ParamAccess.list);
            pManager.AddPointParameter("Right", "R", "Points that are right of the line from P0 to P1", GH_ParamAccess.list);
            pManager.AddPointParameter("On", "O", "Points that are on the line from P0 to P1", GH_ParamAccess.list);
            pManager.AddNumberParameter("Value", "V", "Returns a positive number if it is Left, a negative number if it is Right and 0 if it is on the line", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d P0 = Point3d.Origin;
            Point3d P1 = Point3d.Origin;
            List<Point3d> Pts = new List<Point3d>();

            if (!DA.GetData(0, ref P0)) return;
            if (!DA.GetData(1, ref P1)) return;
            if ((Pts != null) && !DA.GetDataList(2, Pts)) return;

            if (Pts.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Point list is empty");
                return;
            }

            List<Point3d> leftPts = new List<Point3d>();
            List<Point3d> rightPts = new List<Point3d>();
            List<Point3d> onPts = new List<Point3d>();
            List<double> isLeft = new List<double>();

            for (int i = 0; i < Pts.Count; i++)
            {
                double isleft= IsLeft(P0, P1, Pts[i]);
                isLeft.Add(isleft);

                if (isleft > 0.0)
                    leftPts.Add(Pts[i]);
                else if (isleft < 0.0)
                    rightPts.Add(Pts[i]);
                else
                    onPts.Add(Pts[i]);
            }

            DA.SetDataList(0, leftPts);
            DA.SetDataList(1, rightPts);
            DA.SetDataList(2, onPts);
            DA.SetDataList(3, isLeft);
        }

        /// <summary>
        /// Determines if a third point is Left, On or Right of an infinite 2D line defined by two points
        /// </summary>
        /// <returns>Returns a positive number if it is Left, a negative number if it is Right and 0 if it is on the line</returns>
        /// <param name="P0">P0 - the first point</param>
        /// <param name="P1">P1 - the second point</param>
        /// <param name="P2">P2 - the point to test</param>
        private static double IsLeft(Point3d P0, Point3d P1, Point3d P2)
        {
            return ((P1.X - P0.X) * (P2.Y - P0.Y) - (P2.X - P0.X) * (P1.Y - P0.Y));
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.ArePointsLeft_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("32a843ec-34f9-4214-b808-32efc13e77e7"); }
        }
    }
}
