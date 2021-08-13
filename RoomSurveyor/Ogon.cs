using System;
using Rhino.Geometry;


namespace RoomSurveyor

{
    public class Ogon : Polyline
    {
        public Ogon()
        {
        }

        CircularLinkedList<Point3d> rOgon = new CircularLinkedList<Point3d>();
        
    }
}
