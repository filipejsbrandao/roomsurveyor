using System;
using System.Collections.Generic;
using Rhino;
using Rhino.Geometry;


namespace RoomSurveyor
{
    public static class ConvexHull
    {

        /// <summary>
        /// This implements the algorthm described by Lee 83 to determine the Convex Hull of a Simple polygon in O(n)
        /// </summary>
        /// <param name="poly">A simple (non self-intersecting) closed polygon</param>
        /// <returns>A convex polygon</returns>
        public static Polyline Lee83(Polyline poly)
        {
            //NOTE: We could work with 2D points but then we would need to implement a Line class or complicate the implementation since Rhino Line class uses Point3d.

            //ensure the starting point is lower Y largest X
            poly = ShiftLeftToLowerY(poly);

            //List<Point3d> cHull = new List<Point3d>();
            CircularLinkedList<Point3d> cHull = new CircularLinkedList<Point3d>();
            LinkedList<Point3d> polygon = new LinkedList<Point3d>(poly.ToArray());

            if (poly.Count == 3)
            {
                return poly;
            }

            //INICIALIZATION
            //Put two points on the stack
            //int n = polygon.Count;
            LinkedListNode<Point3d> startTop = polygon.First.Next; //top is the top vertex of the stack
            LinkedListNode<Point3d> next = startTop.Next; //next is the next vertex to process or the active vertex
            cHull.AddFirst(polygon.First.Value);
            //cHull.AddLast(polygon.First.Value);
            cHull.AddLast(startTop.Value);
            Node<Point3d> topNode = cHull.Tail;
            Line l = new Line(cHull.Head.Value, topNode.Value);

            do
            {
                if (IsLeft(l.From, l.To, next.Value) <= 0)//Case 1
                {
                    Update(ref topNode, ref next, ref l, ref cHull);
                    CloseLobe(ref topNode, ref next, ref l, ref cHull);
                }
                else //Case 2
                {
                    if (IsLeft(topNode.Value, polygon.First.Value, next.Value) >= 0) //Case 2b1
                    {
                        do
                        {
                            next = next.Next;
                        } while (next != polygon.Last && IsLeft(topNode.Value, polygon.First.Value, next.Value) >= 0);
                    }
                    else //Case 2b2
                    {
                        Line l1 = new Line(topNode.Value, next.Value);
                        CloseLobe(ref topNode, ref next, ref l1, ref cHull);
                        l = l1;
                    }
                }
            } while (next != polygon.Last);

            cHull.AddLast(cHull.Head.Value);

            //int size = cHull.Count;
            //cHull.CopyTo(new Point3d[size], 0);
            //Point3d[] chArray = new Point3d[size];
            
            //return new Polyline(chArray);
            return new Polyline(cHull);
        }
        /// <summary>
        /// <summary>
        /// This routine has as input the top vertex u of the stack, the active
        /// vertex v, and the directed line l = uv.
        /// It checks first to see if the next vertex v' lies in the lobe L(u, v),
        /// i.e., region R0. Upon return, v contains a new active vertex and l
        /// a new current line.
        /// </summary>
        /// <param name="u">the vertex at the top of the stack</param>
        /// <param name="v">the active vertex</param>
        /// <param name="l">the directed line from u to v</param>
        /// <param name="cHull">a stack containing the current convex hull</param>
        /// <param name="poly">the original closed polyline</param>
        private static void CloseLobe(ref Node<Point3d> u, ref LinkedListNode<Point3d> v, ref Line l, ref CircularLinkedList<Point3d> cHull)
        {

            LinkedListNode<Point3d> v1 = v.Next;
            //if v' is not strictly to the right of l and v' is to the left of (v, CW(v))"
            if (IsLeft(l.From, l.To, v1.Value) >= 0 && IsLeft(v.Value, v.Previous.Value, v1.Value) > 0)
            {
                do
                {
                    v1 = v1.Next;
                } while (IsLeft(l.From, l.To, v1.Value) >= 0);
                v = v1;
                l = new Line(u.Previous.Value, u.Value);
            }
            else
            {
                cHull.AddLast(v.Value); //adding v to the stack
                u = cHull.Tail; // v is the new top vertex of the stack
                v = v1; // v1 is the new active vertex
            }
        }

        /// <summary>
        /// It updates the stack by deleting vertices from the stack that are not
        /// hull vertices due to the presence of v. Upon return, u contains the
        /// new top vertex of the stack and l the directed line determined by
        /// u and v
        /// </summary>
        /// <returns>  </returns>
        /// <param name="u">the vertex at the top of the stack</param>
        /// <param name="v">the active vertex</param>
        /// <param name="l">the directed line from u to v</param>
        /// <param name="cHull">a stack containing the current convex hull</param>
        private static void Update(ref Node<Point3d> u, ref LinkedListNode<Point3d> v, ref Line l, ref CircularLinkedList<Point3d> cHull)
        {
            Node<Point3d> v1 = u; //Node<Point3d> predV1;

            do 
            {
                if (v1.Equals(cHull.Head))
                    break;
                v1 = v1.Previous;// v1 is the predecessor of u on the stack //This needs to be revised as v1 will always be set to the predecessor of u...
                cHull.RemoveTail(); // Because u was already added to the stack
                if (v1.Equals(cHull.Head))
                    break;
            } while (IsLeft(v1.Previous.Value, v1.Value, v.Value) <= 0);
            u = v1;
            l = new Line(u.Value, v.Value);
        }

        /// <summary>
        /// Shifts the start point of a closed polyline (ie. a polygon) CCW to the point with smaller Y and largest X.
        /// The polygon must be ccw oriented
        /// </summary>
        /// <param name="poly">The polyline to change</param>
        /// <returns>the rotated polyline</returns>
        private static Polyline ShiftLeftToLowerY(Polyline poly)
        {
            int lowerY = 0;
            List<Point3d> points = new List<Point3d>(poly);
            points.RemoveAt(points.Count - 1);
            //points.AddRange(poly);
            double Ycoord = points[0].Y;
            double Xcoord = points[0].X;

            for (int i = 1; i < points.Count; i++)
            {
                if (points[i].Y <= Ycoord)
                {
                    if (points[i].Y < Ycoord)
                    {
                        lowerY = i;
                        Ycoord = points[i].Y;
                        Xcoord = points[i].X;
                    }
                    else if (points[i].Y == Ycoord && points[i].X > Xcoord)
                    {
                        lowerY = i;
                        Ycoord = points[i].Y;
                        Xcoord = points[i].X;
                    }
                }
            }
            poly.Clear();
            poly.AddRange(OrientPolygon.ShiftLeft(points, lowerY));
            poly.Add(poly[0]);

            //shift = Math.Abs(shift) % points.Count;
            return poly;
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
    }
}
