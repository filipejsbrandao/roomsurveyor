
using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace RoomSurveyor

{
    public class RoomSurveyorCategoryIcon : Grasshopper.Kernel.GH_AssemblyPriority
    {
        public override Grasshopper.Kernel.GH_LoadingInstruction PriorityLoad()
        {
            Grasshopper.Instances.ComponentServer.AddCategoryIcon("RoomSurveyor", RoomSurveyor.Properties.Resources.RoomSurveyor_A_Icon);
            Grasshopper.Instances.ComponentServer.AddCategorySymbolName("RoomSurveyor", 'R');
            return Grasshopper.Kernel.GH_LoadingInstruction.Proceed;
        }
    }

    public class RoomSurveyorInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "RoomSurveyor";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                return Properties.Resources.RoomSurveyor_A_Icon;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "This plug-in features two implementations of interactive triangulation algorithms for 2D polygons which can be used to survey rooms. It includes a number of other tools to manipulate 2D polygons";
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
                return "https://filipebrandao.pt";
            }
        }
        public override string AssemblyVersion => "0.7.1.0";
        
        public override string Version => "0.7.1.0";

        public override Bitmap AssemblyIcon => Properties.Resources.RoomSurveyor_A_Icon;

        public override string AssemblyDescription => "Interactive 2D triangulation";

        public override string AssemblyName => "RoomSurveyor";
    }
}
