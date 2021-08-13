﻿using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor

{
    public class IsPointLeft : GH_Component
    {
        public override Grasshopper.Kernel.GH_Exposure Exposure
        {
            get { return GH_Exposure.hidden; }
        }
        public IsPointLeft()
          : base("IsPointLeft", "IsLeft",
            "Determines if the point is left, on or right of the infinite line from P0 to P1 in 2D space",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("First Point", "P0", "The origin point", GH_ParamAccess.item);
            pManager.AddPointParameter("Second Point", "P1", "The second point determines the direction of the line", GH_ParamAccess.item);
            pManager.AddPointParameter("Test Point", "P2", "Point to test", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Relation", "R", "The relative position of the point in relation to the line. If a point is to the right it will return -1 ,if left it will return 1, and 0 if is on the line", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d P0 = Point3d.Origin;
            Point3d P1 = Point3d.Origin;
            Point3d P2 = Point3d.Origin;
            int relation;

            if (!DA.GetData(0, ref P0)) return;
            if (!DA.GetData(1, ref P1)) return;
            if (!DA.GetData(2, ref P2)) return;

            double isLeft = IsLeft(P0, P1, P2);

            if (isLeft < -Rhino.RhinoMath.ZeroTolerance)
            {
                relation = -1;
            }
            else if(isLeft > Rhino.RhinoMath.ZeroTolerance)
            {
                relation = 1;
            }
            else
            {
                relation = 0;
            }

            DA.SetData(0, relation);

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

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.IsPointLeft_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("863c15da-a228-48fc-89fe-5b13f0f02ea2"); }
        }
    }
}