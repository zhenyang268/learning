#include <math.h>

#define PI 3.14159265358979323846f

typedef struct
{
    float ux, uy, uz; // 旋转轴
    float angle;    // 旋转角度
} AxisAngle;

typedef struct
{
    float yaw;      // 偏航角绕Z, 飞机头左右
    float pitch;    // 俯仰角绕Y, 飞机头上下
    float roll;     // 偏转角绕X, 飞机头旋转
} EulerZYX;


class Rotate
{
    /*
        罗德里格任意轴旋转
        R=cosθ⋅I+(1−cosθ)⋅uuT+sinθ⋅[u]x
        ​
        欧拉叫旋转
        R=RXRYRZ
    */
public:
    // 任意轴旋转, 并返回欧拉角用于验证正确性
    EulerZYX RodriguesRotation();
    int EularRotation();

private:
    
};