using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace RoomSurveyorRH6
{
    public class RoomSurveyorRH6Info : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "Room Surveyor";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                return Properties.Resources.RoomSurvey_Icon;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "This plug-in features an implementation of interactive triangulation algorithms for 2D polygons which can be used to survey rooms. It includes a number of other tools to manipulate 2D polygons";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("329f8870-a398-4c20-9bb2-eb97748525ce");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "Filipe JS Brandão";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "http://filipebrandao.pt";
            }
        }
    }
}
