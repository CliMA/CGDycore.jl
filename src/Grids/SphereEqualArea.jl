function EqualAreaSphereGrid()

end

A = new EqualSphArea(P.ToArray(), targetSum, strength);
public class EqualSphArea : GoalObject
{
double targetSum;
public EqualSphArea(Point3d[] pts, double targetSum, double strength)
{
PPos = pts;
Move = new Vector3d[pts.Length];
Weighting = new double[pts.Length];
for(int i = 0;i < pts.Length;i++) Weighting[i] = strength;
this.targetSum = targetSum;
}
public override void Calculate(List < KangarooSolver.Particle > p)
{
//take sets of spherical polygons’ points
//find angles
//get total angle and required change
//for each angle,
//get bisector of the spherical angle at each points
//move angle pt along vector to this
var pts = this.GetCurrentPositions(p);
for(int i = 0;i < pts.Length;i++)Move[i] = Vector3d.Zero;
double[] angle = new double[pts.Length];
double angleSum = 0;
Vector3d v1,v2;
Vector3d[] T1 = new Vector3d[pts.Length];
Vector3d[] T2 = new Vector3d[pts.Length];
for(int i = 0;i < pts.Length;i++)
{
if(i = = 0)
{v1 = pts[i + 1] - pts[i];
v2 = pts[pts.Length - 1] - pts[i];
}
else if(i = = pts.Length - 1)
{v1 = pts[0] - pts[i];
v2 = pts[i - 1] - pts[i];
}
else
{v1 = pts[i + 1] - pts[i];
v2 = pts[i - 1] - pts[i];
}
Vector3d vo = pts[i] - Point3d.Origin;
Vector3d va1 = Vector3d.CrossProduct(v1, vo);
T1[i] = Vector3d.CrossProduct(vo, va1);
Vector3d va2 = Vector3d.CrossProduct(v2, vo);
T2[i] = Vector3d.CrossProduct(vo, va2);
T1[i].Unitize();
T2[i].Unitize();
Plane pl1 = new Plane (pts[i], pts[i] - Point3d.Origin);
angle[i] = Vector3d.VectorAngle(T1[i], T2[i], pl1);
angleSum + = angle[i];
}
double m = targetSum / angleSum;
for(int i = 0;i < pts.Length;i++)
{
Vector3d vb = T1[i] + 1.05 ∗ T2[i];
vb.Unitize();
if(angle[i] > Math.PI)
{Move[i] + = vb ∗ 0.5 ∗ (m - 1.0);
}
else
{Move[i] + = vb ∗ 0.5 ∗ (1.0 - m);
}
}
}
public override object Output(List < KangarooSolver.Particle > p)
{
return null;

