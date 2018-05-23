#ifndef  VECTOR_H
#define VECTOR_H
#include <math.h>
//#include "intpoint.h"

//cstyle
struct Vector3
{
	double X;
	double Y;
	double Z;

	Vector3(double x = 0, double y = 0, double z=0)
	{
		X = x;
		Y = y;
		Z = z;
	}	

	Vector3 operator+(Vector3 value){ return Vector3(X + value.X, Y + value.Y, Z + value.Z); }
	Vector3 operator-(Vector3 value){ return Vector3(X - value.X, Y - value.Y, Z - value.Z); }
	Vector3 operator*(Vector3 value){ return Vector3(X * value.X, Y * value.Y, Z * value.Z); }
	Vector3 operator/(Vector3 value){ return Vector3(X / value.X, Y / value.Y, Z / value.Z); }

	void operator+=(Vector3 value){ X += value.X; Y += value.Y; Z += value.Z; }
	void operator-=(Vector3 value){ X -= value.X; Y -= value.Y; Z -= value.Z; }
	void operator*=(Vector3 value){ X *= value.X; Y *= value.Y; Z *= value.Z; }
	void operator/=(Vector3 value){ X /= value.X; Y /= value.Y; Z /= value.Z; }

	Vector3 operator+(double value){ return Vector3(X + value, Y + value, Z + value); }
	Vector3 operator-(double value){ return Vector3(X - value, Y - value, Z - value); }
	Vector3 operator*(double value){ return Vector3(X * value, Y * value, Z * value); }
	Vector3 operator/(double value){ return Vector3(X / value, Y / value, Z / value); }

	void  operator+=(double value){ X += value; Y += value; Z += value; }
	void  operator-=(double value){ X -= value; Y -= value; Z -= value; }
	void  operator*=(double value){ X *= value; Y *= value; Z *= value; }
	void  operator/=(double value){ X /= value; Y /= value; Z /= value; }

	double Length()
	{
		return sqrt(X*X + Y*Y + Z*Z);
	}

	Vector3 NormalizeSafe()
	{
		Vector3 rst;
		int l = Length();
		if (l == 0)
		{
			rst.X = rst.Y = rst.Z = 0;
		}
		else
		{
			rst.X = (double)X / (double)l;
			rst.Y = (double)Y / (double)l;
			rst.Z = (double)Z / (double)l;
		}
		return rst;
	}

	Vector3 CrossProduct(Vector3 vector)
	{
		return Vector3(Y * vector.Z - Z * vector.Y,
					   Z * vector.X - X * vector.Z,
					   X * vector.Y - Y * vector.X);
	}
	
	float DotProduct(Vector3 vector)
	{
		float rst = (X * vector.X) + (Y * vector.Y) + (Z * vector.Z);
		return (float)rst;
	}

	bool operator==(Vector3 vector)
	{
		if (X != vector.X || Y != vector.Y || Z != vector.Z)
			return false;
		return true;
	}

	bool operator!=(Vector3 vector)
	{
		return !(*this == vector);
	}
};

struct Vector2
{
	double X;
	double Y;

	Vector2(double x = 0, double y = 0)
	{
		X = x;
		Y = y;
	}

	//Vector2(Point p)
	//{
	//	X = (double)p.X;
	//	Y = (double)p.Y;
	//}

	//Point AsPoint()
	//{
	//	return Point((int)X, (int)Y);
	//}

	Vector2 operator+(Vector2 value){ return Vector2(X + value.X, Y + value.Y); }
	Vector2 operator-(Vector2 value){ return Vector2(X - value.X, Y - value.Y); }
	Vector2 operator*(Vector2 value){ return Vector2(X*value.X, Y*value.Y); }
	Vector2 operator/(Vector2 value){ return Vector2(X / value.X, Y / value.Y); }

	void operator+=(Vector2 value){ X += value.X; Y += value.Y; }
	void operator-=(Vector2 value){ X -= value.X; Y -= value.Y; }
	void operator*=(Vector2 value){ X *= value.X; Y *= value.Y; }
	void operator/=(Vector2 value){ X /= value.X; Y /= value.Y; }

	Vector2 operator+(double value){ return Vector2(X + value, Y + value); }
	Vector2 operator-(double value){ return Vector2(X - value, Y - value); }
	Vector2 operator*(double value){ return Vector2(X*value, Y*value); }
	Vector2 operator/(double value){ return Vector2(X / value, Y / value); }

	void  operator+=(double value){ X += value; Y += value; }
	void  operator-=(double value){ X -= value; Y -= value; }
	void  operator*=(double value){ X *= value; Y *= value; }
	void  operator/=(double value){ X /= value; Y /= value; }

	double Length()
	{
		return sqrt(X*X + Y*Y);
	}

	Vector2 NormalizeSafe()
	{
		Vector2 rst;
		int l = Length();
		if (l == 0)
		{
			rst.X = rst.Y = 0;
		}
		else
		{
			rst.X = (double)X / (double)l;
			rst.Y = (double)Y / (double)l;
		}
		return rst;
	}
};

//INLINE Vector2 normalizeSafe(const Point &p)
//{
//	Vector2 rst;
//	int l = vSize(p);
//	if (l == 0)
//	{
//		rst.X = rst.Y = 0;
//	}
//	else
//	{
//		rst.X = (double)p.X / (double)l;
//		rst.Y = (double)p.Y / (double)l;
//	}
//	return rst;
//}

#endif//VECTOR_H