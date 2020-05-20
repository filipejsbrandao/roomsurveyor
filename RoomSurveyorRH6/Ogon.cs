using System;
using Rhino.Geometry;


namespace RoomSurveyorRH6
{
    public class Ogon : Polyline
    {
        public Ogon()
        {
        }

        CircularLinkedList<Point3d> rOgon = new CircularLinkedList<Point3d>();
        
    }
}
