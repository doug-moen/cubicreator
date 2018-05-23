/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#ifndef FLOAT_POINT_H
#define FLOAT_POINT_H

/*
Floating point 3D points are used during model loading as 3D vectors.
They represent millimeters in 3D space.
*/

#include "intpoint.h"

#include <stdint.h>
#include <math.h>

class FPoint3
{
public:
    double x,y,z;
    FPoint3() {}
    FPoint3(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
	//FPoint3(const Point3& p) : x(p.x*.001), y(p.y*.001), z(p.z*.001) {}

    FPoint3 operator+(const FPoint3& p) const { return FPoint3(x+p.x, y+p.y, z+p.z); }
    FPoint3 operator-(const FPoint3& p) const { return FPoint3(x-p.x, y-p.y, z-p.z); }
    FPoint3 operator*(const double f) const { return FPoint3(x*f, y*f, z*f); }
    FPoint3 operator/(const double f) const { return FPoint3(x/f, y/f, z/f); }
    
    FPoint3& operator += (const FPoint3& p) { x += p.x; y += p.y; z += p.z; return *this; }
    FPoint3& operator -= (const FPoint3& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
    
    bool operator==(FPoint3& p) const { return x==p.x&&y==p.y&&z==p.z; }
    bool operator!=(FPoint3& p) const { return x!=p.x||y!=p.y||z!=p.z; }
    
    double max()
    {
        if (x > y && x > z) return x;
        if (y > z) return y;
        return z;
    }
    
    bool testLength(double len)
    {
        return vSize2() <= len*len;
    }
    
    double vSize2()
    {
        return x*x+y*y+z*z;
    }
    
    double vSize()
    {
        return sqrt(vSize2());
    }

	//inline FPoint3 normalized() //const
	//{
	//	return (*this) / vSize();
	//}

	//FPoint3 cross(const FPoint3& p) const
	//{
	//	return FPoint3(
	//		y*p.z - z*p.y,
	//		z*p.x - x*p.z,
	//		x*p.y - y*p.x);
	//}

	//static FPoint3 cross(const Point3& a, const Point3& b)
	//{
	//	return FPoint3(a).cross(FPoint3(b));
		//        FPoint3(
		//            a.y*b.z-a.z*b.y,
		//            a.z*b.x-a.x*b.z,
		//            a.x*b.y-a.y*b.x);
	//}

	//Point3 toPoint3()
	//{
	//	return Point3(x * 1000, y * 1000, z * 1000);
	//}
};

//inline float operator*(FPoint3 lhs, const FPoint3& rhs) 
//{
//	return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
//}

class FPoint2
{
public:
    double x;
	double y;
    
    FPoint2(double _x=0, double _y=0): x(_x), y(_y) {}
    
    FPoint2 operator+(const FPoint2& p) const { return FPoint2(x+p.x, y+p.y); }
    FPoint2 operator-(const FPoint2& p) const { return FPoint2(x-p.x, y-p.y); }
	FPoint2 operator*(const FPoint2& p) const { return FPoint2(x*p.x, y*p.y); }
    FPoint2 operator/(const FPoint2& p) const { return FPoint2(x/p.x, y/p.y); }

	FPoint2 operator+(const double f) const { return FPoint2(x+f, y+f); }
    FPoint2 operator-(const double f) const { return FPoint2(x-f, y-f); }
    FPoint2 operator*(const double f) const { return FPoint2(x*f, y*f); }
    FPoint2 operator/(const double f) const { return FPoint2(x/f, y/f); }
    
    FPoint2& operator += (const FPoint2& p) { x += p.x; y += p.y; return *this; }
    FPoint2& operator -= (const FPoint2& p) { x -= p.x; y -= p.y; return *this; }
	FPoint2& operator *= (const FPoint2& p) { x *= p.x; y *= p.y; return *this; }
	FPoint2& operator /= (const FPoint2& p) { x /= p.x; y /= p.y; return *this; }

	FPoint2& operator += (double value) { x += value; y += value; return *this; }
    FPoint2& operator -= (double value) { x -= value; y -= value; return *this; }
	FPoint2& operator *= (double value) { x *= value; y *= value; return *this; }
	FPoint2& operator /= (double value) { x /= value; y /= value; return *this; }
    
    bool operator==(FPoint2& p) const { return x==p.x&&y==p.y;}
    bool operator!=(FPoint2& p) const { return x!=p.x||y!=p.y; }
    
    double max()
    {
        if (x > y) 
			return x;
		return y;        
    }
    
    bool testLength(double len)
    {
        return vSize2() <= len*len;
    }
    
    double vSize2()
    {
        return x*x+y*y;
    }
    
    double vSize()
    {
        return sqrt(vSize2());
    }

	FPoint2 normalize()
	{
		double size = vSize();
		FPoint2 p(x, y);
		p/=size;
		return p;
	}
};

class FMatrix3x3
{
public:
    double m[3][3];
    
    FMatrix3x3()
    {
        m[0][0] = 1.0;
        m[1][0] = 0.0;
        m[2][0] = 0.0;
        m[0][1] = 0.0;
        m[1][1] = 1.0;
        m[2][1] = 0.0;
        m[0][2] = 0.0;
        m[1][2] = 0.0;
        m[2][2] = 1.0;
    }
    
    Point3 apply(FPoint3 p)// const
    {
        return Point3(
            (p.x * m[0][0] + p.y * m[1][0] + p.z * m[2][0]) * 1000.0,
            (p.x * m[0][1] + p.y * m[1][1] + p.z * m[2][1]) * 1000.0,
            (p.x * m[0][2] + p.y * m[1][2] + p.z * m[2][2]) * 1000.0);
    }
};

#endif//INT_POINT_H
