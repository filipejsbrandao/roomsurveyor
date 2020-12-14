using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Collections;

namespace RoomSurveyor
{
    public class ArePointLeft3D : GH_Component
    {

        public ArePointLeft3D()
          : base("ArePointsLeft", "AreLeft",
            "Determines if the points are left, on or right of the infinite line from P0 to P1 in 3D space. If a plane is not provided the XY plane is used",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("First Point", "P0", "The origin point", GH_ParamAccess.item);
            pManager.AddPointParameter("Second Point", "P1", "The second point determines the direction of the line", GH_ParamAccess.item);
            pManager.AddPointParameter("Points", "Pts", "Points to test", GH_ParamAccess.list);
            pManager.AddPlaneParameter("Plane", "Pl", "Plane the points are in, if no plane is provided the XY plane is used", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Project", "Pr", "Project to plane. Default is false", GH_ParamAccess.item, false);
            pManager[3].Optional = true;
            pManager[4].Optional = true;
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
            Plane plane = Plane.WorldXY;
            List<Point3d> Pts = new List<Point3d>();
            bool project = false;

            if (!DA.GetData(0, ref P0)) return;
            if (!DA.GetData(1, ref P1)) return;
            if ((Pts != null) && !DA.GetDataList(2, Pts)) return;
            DA.GetData(3, ref plane);
            DA.GetData(4, ref project);

            if (Pts.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Point list is empty");
                return;
            }

            if (!plane.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Please provide a valid plane");
                return;
            }

            List<Point3d> leftPts = new List<Point3d>();
            List<Point3d> rightPts = new List<Point3d>();
            List<Point3d> onPts = new List<Point3d>();
            List<double> isLeft = new List<double>();

            List<Point3d> _Pts = new List<Point3d>();
            Point3dList _PtsXY = new Point3dList();
            Point3d _P0 = P0;
            Point3d _P1 = P1;

            Transform transform = Transform.Identity;//We change nothing
            if (!plane.Equals(Plane.WorldXY))//if the provided plane is not WorldXY we need to move the points to WorldXY.
            {
                transform = Transform.PlaneToPlane(plane, Plane.WorldXY);//otherwise we transform the points
            }

            if (project) { 
                _P0 = plane.ClosestPoint(P0); //get the projection of the provided points to the provided plane
                _P1 = plane.ClosestPoint(P1);

                for (int i = 0; i < Pts.Count; i++)
                {
                    _Pts.Add(plane.ClosestPoint(Pts[i]));
                }
                _PtsXY.AddRange(_Pts);
            }
            else
            {
                _PtsXY.AddRange(Pts);
            }
            _P0.Transform(transform);
            _P1.Transform(transform);
            _PtsXY.Transform(transform);

            for (int i = 0; i < _PtsXY.Count; i++)
            {
                double isleft = IsLeft(_P0, _P1, _PtsXY[i]);
                isLeft.Add(isleft);

                if (isleft > Rhino.RhinoMath.ZeroTolerance)
                    if(project)
                        leftPts.Add(_Pts[i]);
                    else
                        leftPts.Add(Pts[i]);
                else if (isleft < -Rhino.RhinoMath.ZeroTolerance)
                    if(project)
                        rightPts.Add(_Pts[i]);
                    else
                        rightPts.Add(Pts[i]);
                else
                    if (project)
                        onPts.Add(_Pts[i]);
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
            get { return new Guid("4be53214-d2dd-4135-8ee3-33da319516e1"); }
        }
    }
}
