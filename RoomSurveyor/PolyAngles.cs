//
// PolyAngles.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2019 ©2019 Filipe Brandao
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
    public class PolyAngles : GH_Component
    {

        public PolyAngles()
          : base("Polygon Angles", "PA",
            "Returns the internal and reflex angles of a polygon at each corner. If the polygon is Counter Clockwise oriented the angles are returned in the CCW order, else the angles are returned in clockwise order",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "The Polygon", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Degrees", "D", "True returns angles in degrees. Default is Radians", GH_ParamAccess.item, false);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Angles", "A", "The internal angles of each corner of the polygon", GH_ParamAccess.list);
            pManager.AddNumberParameter("Reflex Angles", "RA", "The reflex angles of each corner of the polygon", GH_ParamAccess.list);
            pManager.AddPointParameter("Start Point", "P", "The first corner", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            List<double> angles = new List<double>();
            List<double> reflex = new List<double>();
            Point3d startPt = Point3d.Origin;
            bool degrees = false;
            bool isCCW;

            if (!DA.GetData(0, ref curve)) return;
            DA.GetData(1, ref degrees);

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            //To ensure this always works we must move the polygon to the plane XY and then return it to it original position
            curve.TryGetPlane(out Plane curvePlane);

            if (curve.TryGetPolyline(out poly))
            { }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P input expects a polygon to be provided");
                return;
            }

            if (!poly.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Closed polylines must have at least 3 segments");
                return;
            }
            else if (!poly.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }

            startPt = poly[0];
            Transform transform = Transform.PlaneToPlane(curvePlane, Plane.WorldXY);
            isCCW = IsCCW(poly, curvePlane);
            poly.Transform(transform);
            //isCCW = IsCCW(poly, curvePlane);

            for (int i = 0; i < poly.Count - 1; i++)
            {
                if (i == 0)
                {
                    double angle = AngleAtCorner(poly[i], poly[poly.Count - 2], poly[i + 1]);
                    if (!isCCW)
                        angle = 2 * Math.PI - angle;
                    angles.Add(angle);
                    reflex.Add(2 * Math.PI - angle);
                }
                else
                {
                    double angle = AngleAtCorner(poly[i], poly[i - 1], poly[i + 1]);
                    if (!isCCW)
                        angle = 2 * Math.PI - angle;
                    angles.Add(angle);
                    reflex.Add(2 * Math.PI - angle);
                }
                if (degrees) angles[i] = angles[i] * (360.0 / (2 * Math.PI));
                if (degrees) reflex[i] = reflex[i] * (360.0 / (2 * Math.PI));
            }

            DA.SetDataList(0, angles);
            DA.SetDataList(1, reflex);
            DA.SetData(2, startPt);
        }

        /// <summary>
        /// Checks if polygon is CCW oriented.
        /// </summary>
        /// <returns><c>true</c>, if the polygon is Counter Clockwise oriented, <c>false</c> otherwise.</returns>
        /// <param name="poly">A closed polyline that represents the polygon.</param>
        public static bool IsCCW(Polyline poly, Plane plane)
        {
            bool isCCW = false;
            Curve a = poly.ToNurbsCurve();
            if (a.ClosedCurveOrientation(plane) == CurveOrientation.CounterClockwise)
                isCCW = true;

            return isCCW;
        }

        /// <summary>
        /// Calculates the CCW angle at a corner given the previous and the next point
        /// </summary>
        /// <returns>Returns the angle in radians</returns>
        /// <param name="currPt">current corner - the corner</param>
        /// <param name="prevPt">Previous corner - the previous point or the other end of the edge P1 to P0</param>
        /// <param name="nextPt">Next corner - the next point or the other end of the edge P0 tp P2</param>
        private double AngleAtCorner(Point3d currPt, Point3d prevPt, Point3d nextPt)
        {
            double a = currPt.DistanceTo(prevPt);
            double b = currPt.DistanceTo(nextPt);
            double c = prevPt.DistanceTo(nextPt);

            double angle = FindAngleSSS(a, b, c);
            double isLeft = IsLeft(prevPt, nextPt, currPt);

            if (isLeft > 0)
            { // is a concave angle
                angle = 2 * Math.PI - angle;
            }
            else if (Math.Abs(isLeft) < double.Epsilon)
            { // is colinear
                angle = Math.PI;
            }
            return angle;
        }

        /// <summary>
        /// Finds the internal angle of the line with length <paramref name="a"/> with line with length <paramref name="b"/> given the length of the diagonal <paramref name="c"/>.
        /// </summary>
        /// <returns>The internal angle in radians</returns>
        /// <param name="a">The length of first side .</param>
        /// <param name="b">The length of the second side component.</param>
        /// <param name="c">The length of the diagonal</param>
        private double FindAngleSSS(double a, double b, double c)
        {
            return Math.Acos((a * a + b * b - c * c) / (2 * a * b));
        }

        private double IsLeft(Point3d P0, Point3d P1, Point3d P2)
        {
            return ((P1.X - P0.X) * (P2.Y - P0.Y) - (P2.X - P0.X) * (P1.Y - P0.Y));
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.Angles_Icon;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("0d5af829-2a36-4ff2-87e0-6c87be0e6656"); }
        }
    }
}
