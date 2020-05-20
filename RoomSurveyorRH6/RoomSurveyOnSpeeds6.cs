using System;
using System.Collections.Generic;
using System.Linq;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace RoomSurveyorRH6
{
    public class RoomSurveyOnSpeeds6 : GH_Component
    {

        public RoomSurveyOnSpeeds6()
                  : base("RoomSurvey6OnSpeeds", "RS6_OS",
                    "A test component that solves multiple ogon to polygon transformations. This version of algorithm doesn't use Polygon Chain Closure",
                    "RoomSurveyor", "RoomSurvey")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Starting Polygons", "SP", "A list with the user provided poligons", GH_ParamAccess.list);
            pManager.AddCurveParameter("Objective Polygons", "OP", "A list of the polygon we would like to achieve", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Output Polygons", "OP", "A list with the polygons our algorithm is able to produce", GH_ParamAccess.list);
            pManager.AddTextParameter("User Output", "UO", "The user requested diagonals as strings", GH_ParamAccess.tree);
            pManager.AddLineParameter("Diagonals", "D", "A tree with the requested diagonals for each polygon", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Iterations", "I", "The number of diagonals that were requested", GH_ParamAccess.list);
            pManager.AddNumberParameter("Diagonal Lengths", "DL", "The diagonal lengths that were provided to the algorithm", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Here we will provide the algorithm with the starting polygon, catch the outputs and measure the requested diagonals
            DataTree<string> outText = new DataTree<string>();
            DataTree<Line> diagLines = new DataTree<Line>();
            List<Polyline> rebuiltPolies = new List<Polyline>();
            List<int> iterationList = new List<int>();
            DataTree<double> diagLengths = new DataTree<double>();

            List<Curve> ogons = new List<Curve>();
            List<Curve> polies = new List<Curve>();
            if (!DA.GetDataList(0, ogons)) { return; }
            if (!DA.GetDataList(1, polies)) { return; }

            if (!(ogons.Count == polies.Count))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The lists must have the same number of items");
                return;
            }

            for (int i = 0; i < ogons.Count; i++)
            {
                int iterations = 0;
                GH_Path path = new GH_Path(i);
                List<double> lengths = new List<double>();
                List<double> diagonals = new List<double>();

                //------------------------CHECKING THE OGON----------------------------------------------
                if (!ogons[i].IsPlanar())
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The ogon " + i + " must be a planar polyline");
                    return;
                }

                if (ogons[i].TryGetPolyline(out Polyline ogon))
                { }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Couldn't get a polyline from " + i + " polygon of the Starting Polygons list");
                    return;
                }

                if (!ogon.IsValid)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Polyline " + i + " of the Starting Polygons list must have at least 3 segments");
                    return;
                }
                if (!ogon.IsClosed)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Polyline " + i + " of the Starting Polygons list is not closed");
                    return;
                }

                //------------------------CHECKING THE POLY----------------------------------------------
                if (!polies[i].IsPlanar())
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The poly " + i + " must be a planar polyline");
                    return;
                }

                if (polies[i].TryGetPolyline(out Polyline poly))
                { }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Couldn't get a polyline from " + i + " polygon of the Objective Polygons list");
                    return;
                }

                if (!poly.IsValid)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Polyline " + i + " of the Objective Polygons list must have at least 3 segments");
                    return;
                }
                if (!poly.IsClosed)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Polyline " + i + " of the Objective Polygons list is not closed");
                    return;
                }

                //poly = RoomSurvey4.OrientPoly(poly); Not sure if we should do this in theory the DeOgonizer hasn't changed the start point of the polygon or fliped the orientation

                for (int j = 0; j < poly.SegmentCount; j++)
                {
                    lengths.Add(poly.SegmentAt(j).Length);
                }

                bool closed = false;
                do
                {
                    Polyline returnPoly = RoomSurvey(ogon, lengths, diagonals, out List<string> outT, out List<Line> diagL);

                    string[] words = outT[iterations].Split();
                    int from = 0;
                    int to = 0;

                    if (outT[iterations].Substring(0, 21) == "The Polygon is closed")
                    {
                        rebuiltPolies.Add(returnPoly);
                        diagLines.AddRange(diagL, path);
                        diagLengths.AddRange(diagonals, path);
                        outText.AddRange(outT, path);
                        iterationList.Add(iterations);
                        closed = true;
                    }
                    else if (outT[iterations].Substring(0, 40) == "The length of the last provided diagonal")
                    {
                        rebuiltPolies.Add(returnPoly);
                        diagLines.AddRange(diagL, path);
                        diagLengths.AddRange(diagonals, path);
                        outText.AddRange(outT, path);
                        iterationList.Add(iterations);
                        break;
                    }
                    else
                    {
                        try
                        {
                            from = Int32.Parse(words[5]);
                        }
                        catch (FormatException e)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, e.Message);
                        }
                        try
                        {
                            to = Int32.Parse(words[8]);
                        }
                        catch (FormatException e)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, e.Message);
                        }

                        if (RoomSurvey4.IsValidDiagonal(new Line(poly[from], poly[to]), poly))
                        {
                            diagonals.Add(poly[from].DistanceTo(poly[to]));
                        }
                        else
                        {
                            diagonals.Add(-1);
                        }
                    }
                    ++iterations;

                } while (!closed);
            }

            DA.SetDataList(0, rebuiltPolies);
            DA.SetDataTree(1, outText);
            DA.SetDataTree(2, diagLines);
            DA.SetDataList(3, iterationList);
            DA.SetDataTree(4, diagLengths);
        }

        public Polyline RoomSurvey(Polyline poly, List<double> lengths, List<double> diagonals, out List<string> outText, out List<Line> diagLines)
        {
            List<Vector3d> polyVec = new List<Vector3d>();
            outText = new List<string>();
            diagLines = new List<Line>();
            List<int> isTriVec = new List<int>();

            Polyline rebuiltPoly = new Polyline();
            double tol = 0.015; //in case the lenghts are scaled to mm this would became an int
            double error;

            poly = RoomSurvey6.OrientPoly(poly);//Confirm that the polyline is CCW oriented

            //Contruct the vector chain that represents the new polyline
            for (int i = 0; i < poly.Count - 1; i++)
            {
                Vector3d v = new Vector3d(poly[i + 1] - poly[i]);
                v.Unitize();
                v *= lengths[i];
                polyVec.Add(v);
            }

            // verificar se a cadeia de vetores criará um poligono fechado
            if (IsClosed(polyVec, tol))
            {
                // GOT A SOLUTION
                error = ClosingError(polyVec) * 1000;
                rebuiltPoly = RebuildPoly(poly, polyVec);
                outText.Add("The Polygon is closed with a " + error + " mm error");
            }
            else
            {

                //SET THE VECTORS TO NOT FIXED -----------------------------------------------------------------------------------------------------------------------
                //We will use this to lock the corner's angles
                for (int j = 0; j < polyVec.Count; j++)
                    isTriVec.Add(0);

                //SET THE PROCESSING ORDER -------------------------------------------------------------------------------------------------------------------------------------
                //lets create a vector with the order of the points
                List<int> orderedPoints = new List<int>();
                orderedPoints.AddRange(RoomSurvey6.OrderByAngle(poly));

                //Then the function that creates a matrix with all the diagonals of the polygon
                double[,] diagMatrix = RoomSurvey6.UpperMatrix(poly);

                List<int> orderedDiagonals = new List<int>();
                if (IsPolyOrtho(poly))
                {
                    orderedDiagonals.AddRange(RoomSurvey6.DiagonalOrder(diagMatrix));
                }
                else
                {
                    orderedDiagonals.AddRange(RoomSurvey6.DiagonalOrder(diagMatrix, orderedPoints));
                    orderedDiagonals.AddRange(RoomSurvey6.DiagonalOrder(diagMatrix));
                    orderedDiagonals = orderedDiagonals.Distinct().ToList();
                }

                //SUGEST THE FIRST DIAGONALS --------------------------------------------------------------------------------------------------------------------------------------

                for (int c = 0; c < orderedDiagonals.Count; c++)
                {//create a diagonal using the diagonal vector for the user to see and a text intruction indicating the action
                    int k = orderedDiagonals[c];
                    int matrixSize = diagMatrix.GetLength(0);
                    int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                    int j = k + i * (i + 1) / 2 - matrixSize * i;
                    Line diagonal = new Line(poly[i], poly[j]);
                    diagLines.Add(diagonal);
                    outText.Add("Measure the distance from Point " + i + " to Point " + j);
                    if (diagonals.Count > c)
                    {//if the user has supplied a diagonal add it to the LowerMatrix
                        diagMatrix[j, i] = diagonals[c];
                        bool ijIsTri = false;
                        bool jiIsTri = false;
                        /*
                        Vector3d ijChain = new Vector3d();
                        Vector3d jiChain = new Vector3d();
                        //test if the diagonal[c] closes the poligonal chain i to j or the polygonal chain j to i
                        for (int l = i; l < j; l++) { ijChain += polyVec[l]; }
                        for (int n = j; n < polyVec.Count; n++) { jiChain += polyVec[n]; }
                        for (int m = 0; m < i; m++) { jiChain += polyVec[m]; }
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
                        }*/

                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //--------INTERNAL TRIANGLES---------------------------------------------------
                        //the triangles formed by two edges of the polygon and a diagonal
                        //the triangles formed by one edge of the polygon and two diagonals
                        //the internal triangles of the polygon formed by three diagonals
                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if (RoomSurvey6.InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out bool itoj) >= 0 && diagonals[c] > 0)
                        {
                            int P1 = RoomSurvey6.InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out itoj);
                            int P0 = (itoj) ? i : j;
                            int P2 = (itoj) ? j : i;
                            double isleft = IsLeft(poly[P0], poly[P1], poly[P2]);
                            if (!RoomSurvey6.Triangulate(P0, P1, P2, itoj, diagonals[c], isTriVec, polyVec, isleft))
                            {
                                outText.Add("The length of the last provided diagonal " + diagonals[c] + " is longer than the sum of the length of its adjacent walls from " + P0 + " to " + P2);
                                return rebuiltPoly;
                            }
                            if (itoj) { ijIsTri = true; } else { jiIsTri = true; }
                        }

                        //check if the already provided diagonals will triangulate something
                        int storedDiagonal = RoomSurvey6.StoredDiagonal(diagMatrix, isTriVec, out itoj);
                        if (storedDiagonal > 0)
                        {
                            while (storedDiagonal > 0 && RoomSurvey6.Count0(isTriVec.ToArray(), matrixSize - 1, matrixSize, matrixSize) > 1)
                            {
                                i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * storedDiagonal + 0.25));
                                j = storedDiagonal + i * (i + 1) / 2 - matrixSize * i;
                                int Pone = RoomSurvey6.InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out itoj);
                                int Pzero = (itoj) ? i : j;
                                int Ptwo = (itoj) ? j : i;
                                double isL = IsLeft(poly[Pzero], poly[Pone], poly[Ptwo]);
                                if (!RoomSurvey6.Triangulate(Pzero, Pone, Ptwo, itoj, diagMatrix[j, i], isTriVec, polyVec, isL))
                                {
                                    outText.Add("The length of the last provided diagonal " + diagonals[c] + " is longer than the sum of the length of its adjacent walls from " + Pzero + " to " + Ptwo);
                                    return rebuiltPoly;
                                }
                                if (itoj) { ijIsTri = true; } else { jiIsTri = true; }
                                storedDiagonal = RoomSurvey6.StoredDiagonal(diagMatrix, isTriVec, out itoj);
                            }
                        }
                        //Now we need to check if any or both of the polygonal chains became closed with the diagonal
                        if (ijIsTri && !jiIsTri || !ijIsTri && jiIsTri || diagonals[c] < 0 || !ijIsTri && !jiIsTri)
                        {
                            //remove the diagonals between triangulated points from the diagonalOrder list
                            if (ijIsTri && !jiIsTri || !ijIsTri && jiIsTri || diagonals[c] < 0)
                            {
                                RoomSurvey6.RemoveDiagonals(ijIsTri, jiIsTri, i, j, diagMatrix, orderedDiagonals);
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
                            int pattern101 = RoomSurvey6.Sequence101(isTriVec.ToArray(), isTriVec.Count);//HERE ADD CHECKS FOR -1 value
                            if (pattern101 >= 0)
                            {
                                int center = pattern101;
                                List<int> diagonals101 = new List<int>();
                                diagonals101.AddRange(RoomSurvey6.DiagonalOrder101(center, diagMatrix));
                                k = RoomSurvey6.StoredDiagonal(diagMatrix, diagonals101);//returns the k index of the stored diagonal on the lower matrix or -1 if none found
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
                            int P0 = RoomSurvey6.Pattern000(isTriVec.ToArray(), isTriVec.Count, out int P2, out int P1);
                            if (P0 >= 0)
                            {
                                Vector3d p0p1 = new Vector3d(0, 0, 0);
                                Vector3d p1p2 = new Vector3d(0, 0, 0);
                                Vector3d p2p0 = new Vector3d(0, 0, 0);
                                for (int l = P0; l < P1; l++) { p0p1 += polyVec[l]; }
                                for (int n = P1; n < P2; n++) { p1p2 += polyVec[n]; }
                                for (int m = P2; m < polyVec.Count; m++) { p2p0 += polyVec[m]; }
                                for (int p = 0; p < P0; p++) { p2p0 += polyVec[p]; }
                                if (!RoomSurvey6.IsValidTriangle(p0p1.Length, p1p2.Length, p2p0.Length, out double diff))
                                {
                                    goto Skip;
                                }
                                Vector3d newP1P2 = RoomSurvey6.NextVector(Point3d.Origin, p0p1, p0p1.Length, p1p2.Length, p2p0.Length);
                                if (IsLeft(poly[P0], poly[P1], poly[P2]) < 0)
                                {
                                    Vector3d per = p0p1;
                                    per.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                                    newP1P2.Transform(Transform.Mirror(Point3d.Origin, per));
                                }
                                newP1P2.Reverse();
                                double angle1 = RoomSurvey6.AngleBetweenVectors(p1p2, newP1P2, false);
                                for (int n = P1; n < P2; n++)
                                {
                                    Vector3d v = polyVec[n];
                                    v.Transform(Transform.Rotation(angle1, Point3d.Origin));
                                    polyVec[n] = v;
                                }
                                newP1P2.Reverse();
                                p1p2 = newP1P2;
                                Vector3d newP2P0 = RoomSurvey6.NextVector(Point3d.Origin, p1p2, p1p2.Length, p2p0.Length, p0p1.Length);
                                if (IsLeft(poly[P1], poly[P2], poly[P0]) < 0)
                                {
                                    Vector3d per = p1p2;
                                    per.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                                    newP2P0.Transform(Transform.Mirror(Point3d.Origin, per));
                                }
                                newP2P0.Reverse();
                                angle1 = RoomSurvey6.AngleBetweenVectors(p2p0, newP2P0, false);
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
                                outText.Add("The Polygon is closed with a " + error + " mm error");
                                break;
                            }
                        }
                    Skip:
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //OO PATTERN - search for a 00 unique sequence; find the vector from the start to the end of the pattern
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        int startP00 = RoomSurvey6.Pattern00(isTriVec.ToArray(), isTriVec.Count, out int endP00);
                        if (startP00 >= 0)
                        {
                            Vector3d ij = new Vector3d(0, 0, 0);
                            Vector3d ji = new Vector3d(0, 0, 0);
                            for (int n = startP00; n < endP00; n++) { ij += polyVec[n]; }
                            for (int m = endP00; m < polyVec.Count; m++) { ji += polyVec[m]; }
                            for (int l = 0; l < startP00; l++) { ji += polyVec[l]; }
                            double angle = RoomSurvey6.AngleBetweenVectors(ij, ji, false);
                            for (int n = startP00; n < endP00; n++)
                            {
                                Vector3d v = polyVec[n];
                                v.Transform(Transform.Rotation(angle, Point3d.Origin));
                                polyVec[n] = v;
                            }
                            error = ClosingError(polyVec) * 1000;
                            rebuiltPoly = RebuildPoly(poly, polyVec);
                            outText.Add("The Polygon is closed with a " + error + " mm error");
                            break;
                        }
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //If all corners are triangulated, end processing
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if (RoomSurvey6.Count0(isTriVec.ToArray(), matrixSize - 1, matrixSize, matrixSize) <= 1)
                        {
                            error = ClosingError(polyVec) * 1000;
                            rebuiltPoly = RebuildPoly(poly, polyVec);
                            outText.Add("The Polygon is closed with a " + error + " mm error");
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
                        //End Degug method
                        break;
                    }
                }
            }

            return rebuiltPoly;

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
        /// Determines if a third point is Left, On or Right of an infinite 2D line defined by two points
        /// </summary>
        /// <returns>Returns a positive number if it is Left, a negative number if it is Right and 0 if it is on the line</returns>
        /// <param name="P0">P0 - the first point</param>
        /// <param name="P1">P1 - the second point</param>
        /// <param name="P2">P2 - the point to test</param>
        private double IsLeft(Point3d P0, Point3d P1, Point3d P2)
        {
            return ((P1.X - P0.X) * (P2.Y - P0.Y) - (P2.X - P0.X) * (P1.Y - P0.Y));
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.RS_OnSpeeds_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("3dcb1f8b-2d6f-42bc-8170-757845960b43"); }
        }
    }
}
