package main
import (
	"math"
	"math/rand"
)

type sphere struct {
    color [][]float64 // color vector
    r float64
    center [][]float64
    material uint8
    albedo float64 // a degree of how much of the rays hitting the object get absorbed [0, 1)
}

// returns a function that tells if the ray hits this sphere or not
func (sp *sphere) hit(ray, rayOri [][]float64) (float64, float64, bool) {
    oc := matSub(rayOri, sp.center)
    negB := -vecDot(ray, oc)
    bSq := vecDot(ray, ray)
    Dby4 := negB*negB - bSq*(vecDot(oc, oc)-sp.r*sp.r)
    if Dby4 < 0 { // didnt hit
        return  0, 999999999999, false
        }
    t := (negB - math.Sqrt(Dby4))/bSq // no +ve sqrt(D) cuz we want min anyway
    if t < 0 { // ray hitting behind the camera or really close to object ( t < 0.0000001 for really close thing)
        return 0, 999999999999, false
        }
    intersectionPoint := matScalar(ray, t)
    intersectionNormal := matSub(intersectionPoint, sp.center)
    intersectionNormal = vecUnit(intersectionNormal)
    bright := absVal(vecDot(intersectionNormal, vecUnit(vector(-1, -1, -2))))
    return bright, t, true
}

// decides what ray to get. ie reflected or refracted etc based on paameters of the sphere
func (sp *sphere) getRay(ray, point [][]float64) ([][]float64, [][]float64) {
    if sp.material == 0 {
        return sp.diffuse(point)
    } else if sp.material == 1 {
        return sp.reflection(ray, point)
    } else if sp.material == 2 {
        return sp.refraction(ray, point)
    } else {
        return sp.refraction(ray, point)
    }
}

// returns a point on the unit sphere with center 1 unit from the intersection point in the direction of the normal
func (sp *sphere) diffuse(point [][]float64) ([][]float64, [][]float64) {
    unitnormal := vecUnit(matSub(point, sp.center))
    // point = matAdd(point, matScalar(unitnormal, 0.000000001)) // shadow acne ??
    return matAdd(unitnormal, vecUnit(vector(rand.Float64(), rand.Float64(), rand.Float64()))), point
}

// point is the intersection point
func (sp *sphere) reflection(ray, point [][]float64) ([][]float64, [][]float64) {
    normal := matSub(point, sp.center)
    // ray = matScalar(ray, -1)
    rayOri := point
    ray = matAdd(ray, matScalar(normal, -vecDot(ray, vecUnit(normal))*2))
    // ray = matAdd(ray, matScalar(vector(rand.Float64(), rand.Float64(), rand.Float64()), 3)) // fuzzy reflections? didnt work
    return ray, rayOri
}

func (sp *sphere) refraction(ray, point [][]float64) ([][]float64, [][]float64) {
    muglass := 1.3
    muair := 1.0
    radiusvec := matSub(point, sp.center)
    if vecSize(radiusvec) < sp.r {muair, muglass = muglass, muair}
    unitnormal := vecUnit(radiusvec)
    cross := vecCross(ray, unitnormal)
    perplen := -vecSize(cross)*muair/muglass // -ve cuz reasons
    perpdir := vecUnit(vecCross(cross, unitnormal))
    alonglen := -vecDot(ray, unitnormal)
    ray = matAdd(matScalar(unitnormal, -alonglen), matScalar(perpdir, perplen))
    // if vecSize(radiusvec) > sp.r {return sp.refraction(ray, matAdd(matScalar(radiusvec, 0.9999), sp.center))} // shoot ray, not call refraction
    return ray, point // matAdd(matScalar(radiusvec, 1.0001), sp.center)
}
