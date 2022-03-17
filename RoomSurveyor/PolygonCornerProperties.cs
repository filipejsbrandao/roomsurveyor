//
// PolygonCornerProperties.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2020 ©2020 Filipe Brandao
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
    public class PolygonCornerProperties : GH_Component
    {
        public PolygonCornerProperties()
          : base("Polygon Corner Properties", "PCP",
            "Some morphological properties of the polygon that may be used for comparisons between polygons",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed polygon", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Plane", "Pl", "Provide the plane where the polygon lies", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Orientation", "O", "By default the counter-clockwise oriented turns and types are provided. For clockwise turns and type set to false", GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Type", "Ty", "Returns the a string that describes the type of ogon by turns", GH_ParamAccess.item);
            pManager.AddTextParameter("Turn", "Tu", "A string with the turn at each corner of the ogon, L for left turns, R for right turns", GH_ParamAccess.item);
            pManager.AddNumberParameter("Relation", "R", "An integer for each point of the polygon identifying the relation of the point of the previous two points, 1 for a left turn, -1 for a right turn and 0 for a point that is colinear with the previous and the next", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            Plane userPlane = Plane.WorldXY;
            bool orientation = true;

            string type; string turn;
            List<int> relation = new List<int>();

            if (!DA.GetData(0, ref curve)) return;
            if (!DA.GetData(1, ref userPlane)) return;
            if (!DA.GetData(2, ref orientation)) return;

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
            if (!poly.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }

            if (curvePlane.Normal.IsParallelTo(userPlane.Normal, 0.002) == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon is not parallel within tolerance to the provided plane");
                return;
            }

            //Now we must orient the polygon to the user desired orientation.
            poly = OrientPolygon.OrientPoly(poly, userPlane, orientation);

            Transform transform = Transform.PlaneToPlane(userPlane, Plane.WorldXY);
            poly.Transform(transform);
            relation = PolyCornerRelation(poly);// add a new return type of list of integers with 0s

            type = PolyType(relation);
            turn = PolyTurns(relation);

            DA.SetData(0, type);
            DA.SetData(1, turn);
            DA.SetDataList(2, relation);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="poly">A polygon as a polyline</param>
        /// <returns></returns>
        private List<int> PolyCornerRelation(Polyline poly)
        {
            List<int> relation = new List<int>();

            //Turn the polyline into a list of points
            for (int i = 0; i < poly.SegmentCount; i++)
            {
                Point3d P0 = (i == 0) ? poly[poly.SegmentCount - 1] : poly[i - 1];
                Point3d P1 = poly[i];
                Point3d P2 = poly[i + 1];

                double isLeft = IsLeft(P0, P1, P2);

                if (isLeft < 0)
                {
                    relation.Add(-1);
                }
                else if (isLeft > 0)
                {
                    relation.Add(1);
                }
                else
                {
                    relation.Add(0);
                }
            }
            return relation;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="relation">An integer for each point of the polygon identifying the relation of the point of the previous two points, 1 for a left turn, -1 for a right turn and 0 for a point that is colinear with the previous and the next.</param>
        /// <returns></returns>
        private string PolyTurns(List<int> relation)
        {
            string turns;

            // Use StringBuilder for concatenation in tight loops.
            var sb = new System.Text.StringBuilder();
            for (int i = 0; i < relation.Count; i++)
            {
                if (relation[i] > 0)
                {
                    sb.Append("L");
                }
                else if (relation[i] < 0)
                {
                    sb.Append("R");
                }
            }
            turns = sb.ToString();

            return turns;
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="relation">An integer for each point of the polygon identifying the relation of the point of the previous two points, 1 for a left turn, -1 for a right turn and 0 for a point that is colinear with the previous and the next.</param>
        /// <returns></returns>
        private string PolyType(List<int> relation)
        {
            int R = 0;
            int L = 0;
            List<int> cornerCount = new List<int>();
            string type;

            List<string> typeList = new List<string>();

            for (int i = 0; i < relation.Count; i++)
            {
                int next = (i + 1) % relation.Count;
                if (relation[i] == 1)
                    L++;
                else
                    R++;
                if (relation[next] != relation[i])
                {
                    int count = (L != 0) ? L : R;
                    cornerCount.Add(count);
                    string s = (L != 0) ? L.ToString() + "L" : R.ToString() + "R";
                    typeList.Add(s);
                    L = 0;
                    R = 0;
                }
                if (i == relation.Count - 1 && relation[next] == relation[i])
                {
                    if (typeList.Count == 0)//Its convex
                    {
                        int count = (L != 0) ? L : R;
                        cornerCount.Add(count);
                        string s = (L != 0) ? L.ToString() + "L" : R.ToString() + "R";
                        typeList.Add(s);
                    }
                    else
                    {
                        int j = 0;
                        int length = typeList[0].Length - 1;
                        try
                        {
                            j = Int32.Parse(typeList[0].Substring(0, length));
                        }
                        catch (FormatException e)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, e.Message);
                        }
                        string s = (L != 0) ? (L + j).ToString() + "L" : (R + j).ToString() + "R";
                        int count = (L != 0) ? L + j : R + j;

                        typeList.Insert(0, s);
                        typeList.RemoveAt(1);
                        cornerCount.Insert(0, count);
                        cornerCount.RemoveAt(1);
                    }

                }
            }

            //int k = 0;
            int largest = 0;
            for (int j = 0; j < cornerCount.Count; j++)
            {
                if (largest < cornerCount[j])
                {
                    largest = cornerCount[j];
                    //k = j;
                }
            }

            List<int> positions = new List<int>();
            for (int i = 0; i < cornerCount.Count; i++)
            {
                if (cornerCount[i] == largest) positions.Add(i);
            }

            List<KeyValuePair<int, int>> scores = new List<KeyValuePair<int, int>>();

            foreach (int p in positions)
            {
                List<int> tempList = ShiftLeft(cornerCount, p);
                int score = 0;
                for (int i = 0; i < tempList.Count / 2; i++)
                    score += (tempList.Count - i) * tempList[i];
                scores.Add(new KeyValuePair<int, int>(p, score));
            }

            scores.Sort(delegate (KeyValuePair<int, int> firstPair, KeyValuePair<int, int> nextPair)
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            }
              );

            int k = scores[0].Key;

            typeList = ShiftLeft(typeList, k);

            var sb = new System.Text.StringBuilder();
            for (int i = 0; i < typeList.Count; i++)
            {
                sb.Append(typeList[i].ToString());
            }

            type = sb.ToString();

            return type;
        }

        /// <summary>
        /// Shift a list of integers left by a given amount
        /// </summary>
        /// <param name="list">The list to shift</param>
        /// <param name="shiftBy">The amount to shift by as an integer</param>
        /// <returns>The shifted list to the left.</returns>
        private List<int> ShiftLeft(List<int> list, int shiftBy)
        {
            if (list.Count <= shiftBy)
            {
                return list;
            }

            var result = list.GetRange(shiftBy, list.Count - shiftBy);
            result.AddRange(list.GetRange(0, shiftBy));

            return result;
        }

        /// <summary>
        /// Shift a list of string left by a given amount
        /// </summary>
        /// <param name="list">The list to shift</param>
        /// <param name="shiftBy">The amount to shift by as an integer</param>
        /// <returns>The shifted list to the left.</returns>
        private List<string> ShiftLeft(List<string> list, int shiftBy)
        {
            if (list.Count <= shiftBy)
            {
                return list;
            }

            var result = list.GetRange(shiftBy, list.Count - shiftBy);
            result.AddRange(list.GetRange(0, shiftBy));

            return result;
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
                return Properties.Resources.CornerProperties_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("452d3103-a014-4376-bc99-e451a73cbbaa"); }
        }
    }
}
