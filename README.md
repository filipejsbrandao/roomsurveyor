# RoomSurveyor

**RoomSurveyor** implements iterative triangulation algorithms to assist users to survey rooms, either orthogonal or non-orthogonal, and automating the drawing of the as-is survey. RoomSurvey components in particular, can be used to deploy an interactive workflow of space survey for non-expert users that mimics the empiric surveying workflows used by architects. This enables the possibility of developing low-key design interfaces for mass-customizable systems where contextual information is required.

It is the result of the PhD research developed by [Filipe JS Brandão](https://filipebrandao.pt) at ISCTE-Instituto Universitário de Lisboa / [ISTAR-IUL](http://istar.iscte-iul.pt) research center. It includes interactive triangulation algorithms [(da Silva Brandão et al. 2020)](http://link.springer.com/article/10.1007/s00004-020-00491-3) to assist the accurate survey of non-convex rooms and a few other utilities to handle 2D Polygons. Two versions of the RoomSurvey algorithm are offered: the first, described in [(Brandao et al. 2019)](http://papers.cumincad.org/cgi-bin/works/paper/ecaadesigradi2019_473), is more effective in rooms with few non-orthogonal corners; a second one, RoomSurvey Strict, privileges shortest diagonals and is more effective in rooms with no orthogonal corners.

## components

* ArePointsLeft
* IsPointLeft
* IsOrthoPolygon
* IsInsideWn
* IsConvex
* RandomConvexPoly
* RandomOgon
* RoomSurvey
* RoomSurveyStrict
* InternalDiagonals
* PolygonDiagonals
* OrientPolygon
* PolygonAngles
* PolygonCornerProperties
* ShiftStartPoint
* RoomTurtle



## Installation
You can install RoomSurveyor from the Rhino builtin package manager (_PackageManager command) on macOS or Windows.
Alternatively you need to download the gha file from [Food4Rhino](https://www.food4rhino.com/app/roomsurveyor), unblock it (if in windows) and place it on the Components folder.
