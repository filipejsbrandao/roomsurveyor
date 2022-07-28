//
// RoomSurvey5.cs
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
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class RoomSurvey5 : GH_Component
    {

        public RoomSurvey5()
          : base("RoomSurvey", "RS",
            "RoomSurvey - An algorithm to assist the survey of room plans using diagonals measured by the user between room corners." +
                " This version of the algorithm is more efective in rooms with a few non-orthogonal corners.", "RoomSurveyor", "RoomSurvey")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A polyline that is garanteed to be a closed, non-intersecting and oriented polygon. Use the CheckPolygon component", GH_ParamAccess.item);
            pManager.AddNumberParameter("Side Lengths (m)", "SL", "A list containing the length in meters of each side in an anticlockwise sequence", GH_ParamAccess.list);
            pManager.AddNumberParameter("Diagonals (m)", "D", "A list containing the length of the diagonals in meters in the sequence requested by the Out of this component", GH_ParamAccess.list);
            pManager.AddPlaneParameter("Plane", "Pl", "The plane the polygon is on", GH_ParamAccess.item, Plane.WorldXY);
            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Output", "Out", "A list of string containing the sequence of needed diagonals to close the polygon", GH_ParamAccess.list);
            pManager.AddCurveParameter("Room", "R", "The surveyed room", GH_ParamAccess.item);
            pManager.AddLineParameter("Diagonals", "D", "Required Diagonals as lines for representation purposes", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Triangulated", "T", "Returns true if the polygon is triangulated within tolerance", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            Plane userPlane = Plane.WorldXY;
            bool triangulated = false;
            List<double> lengths = new List<double>();
            List<double> diagonals = new List<double>();

            if (!DA.GetData(0, ref curve)) return;
            if ((lengths != null) && !DA.GetDataList(1, lengths)) return;
            DA.GetDataList(2, diagonals);
            DA.GetData(3, ref userPlane);

            foreach (double l in lengths)
            {
                //Here I don't protect against negative numbers to allow the user to supply -1 to skip diagonals
                if (Math.Abs(l) < Double.Epsilon)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "All sides of the room must be larger than 0");
                    return;
                }
            }

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            if (!userPlane.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Please provide a valid plane");
                return;
            }

            var events = Rhino.Geometry.Intersect.Intersection.CurveSelf(curve, 0.001);
            if (events.Count != 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The provided polygon is either self-intersecting or is very thin. Tolerance for self-intersections is 0.001");
                return;
            }

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
            if (poly.SegmentCount != lengths.Count)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "List of side lengths must have the same number of entries as the number of polyline sides");
                return;
            }

            if (!curve.IsInPlane(userPlane))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon is not on the provided plane");
                return;
            }

            List<Vector3d> polyVec = new List<Vector3d>();
            List<string> outText = new List<string>();
            List<Line> diagLines = new List<Line>();
            List<int> isTriVec = new List<int>();

            Polyline rebuiltPoly = new Polyline();
            double tol = 0.015; //in case the lenghts are scaled to mm this would became an int
            double error;
            Transform transform = Transform.Identity;//We change nothing
            Transform reverseTrans = Transform.Identity;//We change nothing
            if (userPlane.Normal.IsParallelTo(Plane.WorldXY.Normal, 0.001) == 0)
            {
                transform = Transform.PlaneToPlane(userPlane, Plane.WorldXY);//otherwise we transform the polygon
                reverseTrans = Transform.PlaneToPlane(Plane.WorldXY, userPlane);//And put it back in place
            }

            poly.Transform(transform);
            poly = Triangulation.OrientPoly(poly, Plane.WorldXY);//Confirm that the polyline is CCW oriented on the XY plane

            //Contruct the vector chain that represents the new polyline
            for (int i = 0; i < poly.Count - 1; i++)
            {
                Vector3d v = new Vector3d(poly[i + 1] - poly[i]);
                v.Unitize();
                v *= Math.Abs(lengths[i]);//To prevent negative dimensions
                polyVec.Add(v);
            }

            // verificar se a cadeia de vetores criará um poligono fechado
            if (IsClosed(polyVec, tol))
            {
                // GOT A SOLUTION
                error = ClosingError(polyVec) * 1000;
                rebuiltPoly = RebuildPoly(poly, polyVec);
                rebuiltPoly.Transform(reverseTrans);
                outText.Add("The Polygon is closed with a " + error + " mm error");
                triangulated = true;
            }
            else
            {
                //SET THE VECTORS TO NOT FIXED -----------------------------------------------------------------------------------------------------------------------
                //We will use this to lock the corner's angles
                for (int j = 0; j < polyVec.Count; j++)
                    isTriVec.Add(0);

                //lets create a vector with the order of the points
                List<int> orderedPoints = new List<int>();
                orderedPoints.AddRange(Triangulation.OrderByAngle(poly, false));

                //Then the function that creates a matrix with all the diagonals of the polygon
                double[,] diagMatrix = Triangulation.UpperMatrix(poly);

                //SET THE DIAGONAL PROCESSING ORDER
                List<int> orderedDiagonals = new List<int>();
                if (IsPolyOrtho(poly))
                {
                    orderedDiagonals.AddRange(Triangulation.DiagonalOrder(diagMatrix));
                }
                else
                {
                    orderedDiagonals.AddRange(Triangulation.DiagonalOrder(diagMatrix, orderedPoints));
                    orderedDiagonals.AddRange(Triangulation.DiagonalOrder(diagMatrix));
                    orderedDiagonals = orderedDiagonals.Distinct().ToList();
                }

                //SUGEST THE FIRST DIAGONALS --------------------------------------------------------------------------------------------------------------------------------------
                //The sugestor will propose a diagonal the user can measure
                //If all internal angles are ortho then select the first corner with the longest diagonal
                //If the user cannot measure a diagonal (i.e. its value is -1) the algorithm will cycle to the next

                for (int c = 0; c < orderedDiagonals.Count; c++)
                {//create a diagonal using the diagonal vector for the user to see and a text intruction indicating the action
                    int k = orderedDiagonals[c];
                    int matrixSize = diagMatrix.GetLength(0);
                    int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                    int j = k + i * (i + 1) / 2 - matrixSize * i;
                    Line diagonal = new Line(poly[i], poly[j]);
                    diagonal.Transform(reverseTrans);
                    diagLines.Add(diagonal);
                    outText.Add("Measure the distance from Point " + i + " to Point " + j);
                    if (diagonals.Count > c)
                    {//if the user has supplied a diagonal add it to the LowerMatrix
                        diagMatrix[j, i] = diagonals[c];
                        bool ijIsTri = false;
                        bool jiIsTri = false;
                        Vector3d ijChain = new Vector3d();
                        Vector3d jiChain = new Vector3d();

                        for (int l = i; l < j; l++) { ijChain += polyVec[l]; }
                        for (int n = j; n < polyVec.Count; n++) { jiChain += polyVec[n]; }
                        for (int m = 0; m < i; m++) { jiChain += polyVec[m]; }
                        // A tolerance of 7,5mm causes deviations of around 2.2 cm in a 6.175m diagonal
                        // most recent Laser Distance Meters are rated with a +/- 1.5mm accuracy, in our case user introduced errors are likely higher

                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //--------INTERNAL TRIANGLES---------------------------------------------------
                        //the triangles formed by two edges of the polygon and a diagonal
                        //the triangles formed by one edge of the polygon and two diagonals
                        //the internal triangles of the polygon formed by three diagonals
                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if (Triangulation.InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out bool itoj) >= 0 && diagonals[c] > 0)
                        {
                            int P1 = Triangulation.InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out itoj);
                            int P0 = (itoj) ? i : j;
                            int P2 = (itoj) ? j : i;
                            double isleft = IsLeft(poly[P0], poly[P1], poly[P2]);
                            if (!Triangulation.Triangulate(P0, P1, P2, itoj, diagonals[c], isTriVec, polyVec, isleft))
                            {
                                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The length of the last provided diagonal " + diagonals[c] + " is longer than the sum of the length of its adjacent walls from " + P0 + " to " + P2);
                                outText.Add("The length of the last provided diagonal " + diagonals[c] + " is longer than the sum of the length of its adjacent walls from " + P0 + " to " + P2);
                                DA.SetDataList(0, outText);
                                DA.SetDataList(2, diagLines);
                                return;
                            }
                            if (itoj) { ijIsTri = true; } else { jiIsTri = true; }
                        }

                        //check if the already provided diagonals will triangulate something
                        int storedDiagonal = Triangulation.StoredDiagonal(diagMatrix, isTriVec, out itoj);
                        if (storedDiagonal > 0)
                        {
                            while (storedDiagonal > 0 && Triangulation.Count0(isTriVec.ToArray(), matrixSize - 1, matrixSize, matrixSize) > 1)
                            {
                                i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * storedDiagonal + 0.25));
                                j = storedDiagonal + i * (i + 1) / 2 - matrixSize * i;
                                int Pone = Triangulation.InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out itoj);
                                int Pzero = (itoj) ? i : j;
                                int Ptwo = (itoj) ? j : i;
                                double isL = IsLeft(poly[Pzero], poly[Pone], poly[Ptwo]);
                                if (!Triangulation.Triangulate(Pzero, Pone, Ptwo, itoj, diagMatrix[j, i], isTriVec, polyVec, isL))
                                {
                                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The length of the last provided diagonal " + diagMatrix[j, i] + " is longer than the sum of the length of its adjacent walls from " + Pzero + " to " + Ptwo);
                                    outText.Add("The length of the last provided diagonal " + diagMatrix[j, i] + " is longer than the sum of the length of its adjacent walls from " + Pzero + " to " + Ptwo);
                                    DA.SetDataList(0, outText);
                                    DA.SetDataList(2, diagLines);
                                    return;
                                }
                                if (itoj) { ijIsTri = true; } else { jiIsTri = true; }
                                storedDiagonal = Triangulation.StoredDiagonal(diagMatrix, isTriVec, out itoj);
                            }
                        }
                        //test if the diagonal[c] closes the poligonal chain i to j or the polygonal chain j to i
                        //Now we need to check if any or both of the polygonal chains became closed with the diagonal
                        if (Math.Abs(ijChain.Length - diagonals[c]) < tol / 2 && diagonals[c] > 0)
                        {
                            ijIsTri = true;//lock the vectors from i+1 to j
                            for (int n = i + 1; n < j; n++)
                            {
                                isTriVec[n] = 1;
                            }
                        }

                        if (Math.Abs(jiChain.Length - diagonals[c]) < tol / 2 && diagonals[c] > 0)
                        {
                            jiIsTri = true;//lock the vectors of the polygonal chain from j+1 to i
                            for (int n = j + 1; n < polyVec.Count; n++)
                            {
                                isTriVec[n] = 1;
                            }

                            for (int m = 0; m < i; m++)
                            {
                                isTriVec[m] = 1;
                            }
                        }

                        if (ijIsTri && !jiIsTri || !ijIsTri && jiIsTri || diagonals[c] < 0 || !ijIsTri && !jiIsTri)
                        {
                            //remove the diagonals between triangulated points from the diagonalOrder list
                            if (ijIsTri && !jiIsTri || !ijIsTri && jiIsTri || diagonals[c] < 0)
                            {
                                Triangulation.RemoveDiagonals(ijIsTri, jiIsTri, i, j, diagMatrix, orderedDiagonals);
                                if (orderedDiagonals.Count <= outText.Count && !triangulated)
                                {
                                    outText.Add("No more diagonals to request");
                                }
                            }
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            //--------101 PATTERN-------------
                            //101 paterns can speedup the triangulation and help to triangulate concave corners. The aim is always to fix the 0 in the sequence.
                            //we first create an ordered list of valid diagonals from large to small. 
                            //Request the user the largest valid diagonal and if not possible then request one of the remaining ones.
                            //if no valid diagonals exist then continue.
                            //There can be multiple 101 patterns, the key is to find the first zero on both sides: 0...101.0
                            //if the center 0 is i we need to find the i-nth 0 and the i+nth 0 and all the possible diagonals inside.
                            //order the diagonals by size
                            //CURRENTLY DIAGONALORDER101 ONLY RETURNS 4 DIAGONALS
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            int pattern101 = Triangulation.Sequence101(isTriVec.ToArray(), isTriVec.Count);//HERE ADD CHECKS FOR -1 value
                            if (pattern101 >= 0)
                            {
                                int center = pattern101;
                                List<int> diagonals101 = new List<int>();
                                diagonals101.AddRange(Triangulation.DiagonalOrder101(center, diagMatrix));
                                k = Triangulation.StoredDiagonal(diagMatrix, diagonals101);//returns the k index of the stored diagonal on the lower matrix or -1 if none found
                                if (k < 0)
                                {//determine if the first diagonal is on the list (but not before the current c) and remove it
                                    if (diagonals101.Any())
                                    {
                                        int found = orderedDiagonals.IndexOf(diagonals101[0], c + 1);
                                        if (found != -1)
                                        {
                                            orderedDiagonals.RemoveAt(found);
                                        }//insert it on the next c
                                        orderedDiagonals.Insert(c + 1, diagonals101[0]);
                                    }
                                }
                            }
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            //OOO Pattern - search for a 000 unique sequence; find the two intersections of the last not fixed edges and select the correct one
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            int P0 = Triangulation.Pattern000(isTriVec.ToArray(), isTriVec.Count, out int P2, out int P1);
                            if (P0 >= 0)
                            {
                                Vector3d p0p1 = new Vector3d(0, 0, 0);
                                Vector3d p1p2 = new Vector3d(0, 0, 0);
                                Vector3d p2p0 = new Vector3d(0, 0, 0);
                                for (int l = P0; l < P1; l++) { p0p1 += polyVec[l]; }
                                for (int n = P1; n < P2; n++) { p1p2 += polyVec[n]; }
                                for (int m = P2; m < polyVec.Count; m++) { p2p0 += polyVec[m]; }
                                for (int p = 0; p < P0; p++) { p2p0 += polyVec[p]; }
                                if (!Triangulation.IsValidTriangle(p0p1.Length, p1p2.Length, p2p0.Length, out double diff))
                                {
                                    goto Skip;
                                }
                                Vector3d newP1P2 = Triangulation.NextVector(Point3d.Origin, p0p1, p0p1.Length, p1p2.Length, p2p0.Length);
                                if (IsLeft(poly[P0], poly[P1], poly[P2]) < 0)
                                {
                                    Vector3d per = p0p1;
                                    per.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                                    newP1P2.Transform(Transform.Mirror(Point3d.Origin, per));
                                }
                                newP1P2.Reverse();
                                double angle1 = Triangulation.AngleBetweenVectors(p1p2, newP1P2, false);
                                for (int n = P1; n < P2; n++)
                                {
                                    Vector3d v = polyVec[n];
                                    v.Transform(Transform.Rotation(angle1, Point3d.Origin));
                                    polyVec[n] = v;
                                }
                                newP1P2.Reverse();
                                p1p2 = newP1P2;
                                Vector3d newP2P0 = Triangulation.NextVector(Point3d.Origin, p1p2, p1p2.Length, p2p0.Length, p0p1.Length);
                                if (IsLeft(poly[P1], poly[P2], poly[P0]) < 0)
                                {
                                    Vector3d per = p1p2;
                                    per.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                                    newP2P0.Transform(Transform.Mirror(Point3d.Origin, per));
                                }
                                newP2P0.Reverse();
                                angle1 = Triangulation.AngleBetweenVectors(p2p0, newP2P0, false);
                                for (int n = P2; n < polyVec.Count; n++)
                                {
                                    Vector3d v = polyVec[n];
                                    v.Transform(Transform.Rotation(angle1, Point3d.Origin));
                                    polyVec[n] = v;
                                }
                                for (int m = 0; m < P0; m++)
                                {
                                    Vector3d v = polyVec[m];
                                    v.Transform(Transform.Rotation(angle1, Point3d.Origin));
                                    polyVec[m] = v;
                                }
                                isTriVec[P0] = 1;
                                isTriVec[P1] = 1;
                                isTriVec[P2] = 1;
                                error = ClosingError(polyVec) * 1000;
                                rebuiltPoly = RebuildPoly(poly, polyVec);
                                rebuiltPoly.Transform(reverseTrans);
                                outText.Add("The Polygon is closed with a " + error + " mm error");
                                triangulated = true;
                                break;
                            }
                        }
                    Skip:
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //OO PATTERN - search for a 00 unique sequence; find the vector from the start to the end of the pattern
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        int startP00 = Triangulation.Pattern00(isTriVec.ToArray(), isTriVec.Count, out int endP00);
                        if (startP00 >= 0)
                        {
                            Vector3d ij = new Vector3d(0, 0, 0);
                            Vector3d ji = new Vector3d(0, 0, 0);
                            for (int n = startP00; n < endP00; n++) { ij += polyVec[n]; }
                            for (int m = endP00; m < polyVec.Count; m++) { ji += polyVec[m]; }
                            for (int l = 0; l < startP00; l++) { ji += polyVec[l]; }
                            double angle = Triangulation.AngleBetweenVectors(ij, ji, false);
                            for (int n = startP00; n < endP00; n++)
                            {
                                Vector3d v = polyVec[n];
                                v.Transform(Transform.Rotation(angle, Point3d.Origin));
                                polyVec[n] = v;
                            }
                            error = ClosingError(polyVec) * 1000;
                            rebuiltPoly = RebuildPoly(poly, polyVec);
                            rebuiltPoly.Transform(reverseTrans);
                            outText.Add("The Polygon is closed with a " + error + " mm error");
                            triangulated = true;
                            break;
                        }
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //If all corners are triangulated, end processing
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if (Triangulation.Count0(isTriVec.ToArray(), matrixSize - 1, matrixSize, matrixSize) <= 1)
                        {
                            error = ClosingError(polyVec) * 1000;
                            rebuiltPoly = RebuildPoly(poly, polyVec);
                            rebuiltPoly.Transform(reverseTrans);
                            outText.Add("The Polygon is closed with a " + error + " mm error");
                            triangulated = true;
                            break;
                        }
                    }
                    else
                    {
                        //Debug Method
                        rebuiltPoly.Add(poly[0]);
                        Point3d pt = new Point3d();

                        for (int m = 0; m < polyVec.Count; m++)
                        {
                            pt = rebuiltPoly[m] + polyVec[m];
                            rebuiltPoly.Add(pt);
                        }
                        rebuiltPoly.Transform(reverseTrans);
                        //End Degug method
                        break;
                    }
                }
            }
            DA.SetDataList(0, outText);
            DA.SetData(1, rebuiltPoly);
            DA.SetDataList(2, diagLines);
            DA.SetData(3, triangulated);

        }

        /// <summary>
        /// Returns the tolerance required for closing the polygon
        /// </summary>
        /// <returns>The effective tolerance in mm.</returns>
        /// <param name="vectors">A List of vectors that create the polygon.</param>
        private double ClosingError(List<Vector3d> vectors)
        {
            Vector3d s = new Vector3d();
            foreach (Vector3d v in vectors) s += v;
            double error = s.Length;
            return error;
        }

        /// <summary>
        /// Rebuilds a polyline with a list of <paramref name="vectors"/>
        /// </summary>
        /// <returns>The poly.</returns>
        /// <param name="polyline">Polyline.</param>
        /// <param name="vectors">Vectors.</param>
        private Polyline RebuildPoly(Polyline polyline, List<Vector3d> vectors)
        {
            Polyline newPoly = new Polyline();
            newPoly.Add(polyline[0]);
            Point3d p = new Point3d();
            for (int j = 0; j < vectors.Count; j++)
            {
                if (j == vectors.Count - 1)
                {
                    p = newPoly[0];
                    newPoly.Add(p);
                }
                else
                {
                    p = newPoly[j] + vectors[j];
                    newPoly.Add(p);
                }
            }
            return newPoly;
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
        /// Checks if a list of vectors creates a closed polygon within a given tolerance.
        /// </summary>
        /// <returns><c>true</c>, if the resultant vector is smaller than the given <paramref name="tolerance"/>, <c>false</c> otherwise.</returns>
        /// <param name="vectors">A list contaning the vectors.</param>
        /// <param name="tolerance">Tolerance.</param>
        private bool IsClosed(List<Vector3d> vectors, double tolerance)
        {
            Vector3d s = new Vector3d();
            foreach (Vector3d v in vectors) s += v;
            //IsTiny method takes a double
            if (s.IsTiny(tolerance) == true || s.IsZero)
                return true;
            else
                return false;

        }
        /// <summary>
        /// Checks if all the angles in polyline are orthogonal, either 90 or 270.
        /// </summary>
        /// <returns><c>true</c>, if an ortho-polygon was input, <c>false</c> otherwise.</returns>
        /// <param name="polyline">Polyline.</param>
        private bool IsPolyOrtho(Polyline polyline)
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

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.RoomSurvey_Icon;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("5d5cd560-e468-42f1-bfad-1e36b3203427"); }
        }
    }
}
