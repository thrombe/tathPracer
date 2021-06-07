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
func (sp *sphere) hit(ray, rayOri [][]float64) (float64, bool) {
    oc := matSub(rayOri, sp.center)
    negB := -vecDot(ray, oc)
    bSq := vecDot(ray, ray)
    Dby4 := negB*negB - bSq*(vecDot(oc, oc)-sp.r*sp.r)
    if Dby4 < 0 { // didnt hit
        return  999999999999, false
        }
    t := (negB - math.Sqrt(Dby4))/bSq // no +ve sqrt(D) cuz we want min anyway
    if t < 0 { // ray hitting behind the camera or really close to object ( t < 0.0000001 for really close thing)
        return 999999999999, false
        }
    return t, true
}

// decides what ray to get. ie reflected or refracted etc based on paameters of the sphere
func (sp *sphere) getRay(ray, rayOri [][]float64, t float64) ([][]float64, [][]float64) {
    if sp.material == 0 {
        rayOri = matAdd(rayOri, matScalar(ray, t))
        return sp.diffuse(rayOri)
    } else if sp.material == 1 {
        rayOri = matAdd(rayOri, matScalar(ray, t))
        return sp.reflection(ray, rayOri)
    } else if sp.material == 2 {
        return sp.refraction(ray, rayOri, t)
    } else {
        return sp.refraction(ray, rayOri, t)
    }
}

// returns a point on the unit sphere with center 1 unit from the intersection point in the direction of the normal
func (sp *sphere) diffuse(rayOri [][]float64) ([][]float64, [][]float64) {
    unitnormal := vecUnit(matSub(rayOri, sp.center))
    return matAdd(unitnormal, vecUnit(vector(rand.Float64(), rand.Float64(), rand.Float64()))), rayOri
}

// point is the intersection point
func (sp *sphere) reflection(ray, rayOri [][]float64) ([][]float64, [][]float64) {
    normal := matSub(rayOri, sp.center)
    ray = matAdd(ray, matScalar(normal, -vecDot(ray, vecUnit(normal))*2))
    // ray = matAdd(ray, matScalar(vector(rand.Float64(), rand.Float64(), rand.Float64()),  vecSize(ray)/5)) // fuzzy reflections? didnt work
    return ray, rayOri
}

// BROKEN
func (sp *sphere) refraction(ray, rayOri [][]float64, t float64) ([][]float64, [][]float64) {
    muglass := 1.5
    muair := 1.0
    rayOri = matAdd(rayOri, matScalar(ray, t*1.00001))
    radiusvec := matSub(rayOri, sp.center)
    // if vecDot(radiusvec, radiusvec) < sp.r*sp.r {muair, muglass, radiusvec = muglass, muair, matScalar(radiusvec, -1)}
    unitnormal := vecUnit(radiusvec)
    raypall := matScalar(unitnormal, vecDot(ray, unitnormal))
    pccrossn := vecCross(matScalar(ray, -1), unitnormal)
    oriperp := matScalar(vecUnit(vecCross(pccrossn, unitnormal)), vecSize(pccrossn))
    rayperp := matScalar(oriperp, muair/muglass)
    ray = matAdd(raypall, rayperp)

    t, hit := sp.hit(ray, rayOri)
    if !hit {return ray, rayOri}
    rayOri = matAdd(rayOri, matScalar(ray, t*1.00001))
    radiusvec = matSub(sp.center, rayOri)
    unitnormal = vecUnit(radiusvec)
    raypall = matScalar(unitnormal, vecDot(ray, unitnormal))
    pccrossn = vecCross(matScalar(ray, -1), unitnormal)
    oriperp = matScalar(vecUnit(vecCross(pccrossn, unitnormal)), vecSize(pccrossn))
    rayperp = matScalar(oriperp, muglass/muair)
    ray = matAdd(raypall, rayperp)
    return ray, rayOri
}

func (sp *sphere) getCol(color [][]float64) [][]float64 {
    if sp.material == 0 {return sp.color}
    if sp.material == 1 {return color}
    if sp.material == 2 {return color}
    return color
}