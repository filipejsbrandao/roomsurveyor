//
// IsPointLeft3D.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2021 ©2021 Filipe Brandao
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class IsPointLeft3D : GH_Component
    {

        public IsPointLeft3D()
          : base("IsPointLeft3D", "IsLeft",
            "Determines if the point is left, on or right of the infinite line from P0 to P1 in 3D space. If a plane is not provided the XY plane is used",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("First Point", "P0", "The origin point", GH_ParamAccess.item);
            pManager.AddPointParameter("Second Point", "P1", "The second point determines the direction of the line", GH_ParamAccess.item);
            pManager.AddPointParameter("Test Point", "P2", "Point to test", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Plane", "Pl", "Plane the points are in, if no plane is provided the XY plane is used", GH_ParamAccess.item);
            pManager[3].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Relation", "R", "The relative position of the point in relation to the line. (right = -1 ,left = 1, coincident = 0)", GH_ParamAccess.item);
            pManager.AddPointParameter("First point'", "P0'", "The first point projected to the plane", GH_ParamAccess.item);
            pManager.AddPointParameter("Second point'", "P1'", "The first point projected to the plane", GH_ParamAccess.item);
            pManager.AddPointParameter("Test point'", "P2'", "The first point projected to the plane", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d P0 = Point3d.Origin;
            Point3d P1 = Point3d.Origin;
            Point3d P2 = Point3d.Origin;
            Plane plane = Plane.WorldXY;

            int relation;

            if (!DA.GetData(0, ref P0)) return;
            if (!DA.GetData(1, ref P1)) return;
            if (!DA.GetData(2, ref P2)) return;
            DA.GetData(3, ref plane);

            if (!plane.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Please provide a valid plane");
                return;
            }

            if (!P0.IsValid || !P1.IsValid || !P2.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Please provide valid points");
                return;
            }

            Point3d _P0 = plane.ClosestPoint(P0); //get the projection of the provided points to the provided plane
            Point3d _P1 = plane.ClosestPoint(P1);
            Point3d _P2 = plane.ClosestPoint(P2);

            DA.SetData(1, _P0);
            DA.SetData(2, _P1);
            DA.SetData(3, _P2);

            Transform transform = Transform.Identity;//We change nothing
            if (!plane.Equals(Plane.WorldXY))//if the provided plane is not WorldXY we need to move the points to WorldXY.
            {
                transform = Transform.PlaneToPlane(plane, Plane.WorldXY);//otherwise we transform the points
            }

            _P0.Transform(transform);
            _P1.Transform(transform);
            _P2.Transform(transform);

            double isLeft = IsLeft(_P0, _P1, _P2);

            if (isLeft < - Rhino.RhinoMath.ZeroTolerance)
            {
                relation = -1;
            }
            else if (isLeft > Rhino.RhinoMath.ZeroTolerance)
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
            get { return new Guid("b8a44d9d-4251-48b7-b169-d9e258d5886a"); }
        }
    }
}
