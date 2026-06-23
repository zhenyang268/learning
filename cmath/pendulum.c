#include <stdio.h>
#include <math.h>
#include <raylib.h>

#define WIDTH 800
#define HEIGHT 600
#define G 9.81f
#define L1 120.f
#define L2 120.f
#define M1 1.f
#define M2 1.f
#define TRAIL_LENGTH 800     // 轨迹点数
#define SUBSTEPS 4            // 每帧物理子步数

typedef struct Pendulum
{
    float th1; // 摆1的角度（从竖直向下=0，逆时针为正）
    float th2; // 摆2的角度
    float w1;  // 摆1的角速度
    float w2;  // 摆2的角速度
} Pendulum;

// 双摆运动方程的导数
// 参考标准 Lagrangian 推导
void Derivs(Pendulum p, float *dth1, float *dth2, float *dw1, float *dw2)
{
    float sin1 = sinf(p.th1), cos1 = cosf(p.th1);
    float sin2 = sinf(p.th2), cos2 = cosf(p.th2);
    float dth = p.th1 - p.th2;
    float sinD = sinf(dth), cosD = cosf(dth);

    float den = M1 + M2 * sinD * sinD;  // 公共分母因子

    // θ₁'' 的分子
    float num1 = M2 * G * sin2 * cosD
               - M2 * sinD * (L1 * p.w1 * p.w1 * cosD + L2 * p.w2 * p.w2)
               - (M1 + M2) * G * sin1;

    // θ₂'' 的分子
    float num2 = (M1 + M2) * (L1 * p.w1 * p.w1 * sinD - G * sin2 + G * sin1 * cosD)
               + M2 * L2 * p.w2 * p.w2 * sinD * cosD;

    *dth1 = p.w1;
    *dth2 = p.w2;
    *dw1 = num1 / (L1 * den);
    *dw2 = num2 / (L2 * den);
}

// 四阶龙格库塔法（RK4）单步积分
void Step(Pendulum *p, float dt)
{
    float k1t1, k1t2, k1w1, k1w2;
    float k2t1, k2t2, k2w1, k2w2;
    float k3t1, k3t2, k3w1, k3w2;
    float k4t1, k4t2, k4w1, k4w2;
    Pendulum tmp;

    // k1
    Derivs(*p, &k1t1, &k1t2, &k1w1, &k1w2);

    // k2
    tmp.th1 = p->th1 + k1t1 * dt * 0.5f;
    tmp.th2 = p->th2 + k1t2 * dt * 0.5f;
    tmp.w1  = p->w1  + k1w1 * dt * 0.5f;
    tmp.w2  = p->w2  + k1w2 * dt * 0.5f;
    Derivs(tmp, &k2t1, &k2t2, &k2w1, &k2w2);

    // k3
    tmp.th1 = p->th1 + k2t1 * dt * 0.5f;
    tmp.th2 = p->th2 + k2t2 * dt * 0.5f;
    tmp.w1  = p->w1  + k2w1 * dt * 0.5f;
    tmp.w2  = p->w2  + k2w2 * dt * 0.5f;
    Derivs(tmp, &k3t1, &k3t2, &k3w1, &k3w2);

    // k4
    tmp.th1 = p->th1 + k3t1 * dt;
    tmp.th2 = p->th2 + k3t2 * dt;
    tmp.w1  = p->w1  + k3w1 * dt;
    tmp.w2  = p->w2  + k3w2 * dt;
    Derivs(tmp, &k4t1, &k4t2, &k4w1, &k4w2);

    // 加权求和
    p->th1 += dt / 6.0f * (k1t1 + 2.0f * k2t1 + 2.0f * k3t1 + k4t1);
    p->th2 += dt / 6.0f * (k1t2 + 2.0f * k2t2 + 2.0f * k3t2 + k4t2);
    p->w1  += dt / 6.0f * (k1w1 + 2.0f * k2w1 + 2.0f * k3w1 + k4w1);
    p->w2  += dt / 6.0f * (k1w2 + 2.0f * k2w2 + 2.0f * k3w2 + k4w2);
}

// 计算摆杆末端坐标
Vector2 GetBobPos(float th, Vector2 origin, float length)
{
    return (Vector2){origin.x + length * sinf(th), origin.y + length * cosf(th)};
}

// 绘制摆
void DrawPendulum(Vector2 origin, float th1, float th2)
{
    Vector2 pos1 = GetBobPos(th1, origin, L1);
    Vector2 pos2 = GetBobPos(th2, pos1, L2);

    // 摆杆
    DrawLineEx(origin, pos1, 3.0f, DARKGRAY);
    DrawLineEx(pos1, pos2, 3.0f, DARKGRAY);

    // 固定点
    DrawCircleV(origin, 6, DARKGRAY);

    // 摆球
    DrawCircleV(pos1, 12, RED);
    DrawCircleV(pos2, 12, BLUE);
}

// 绘制轨迹
void DrawTrail(Vector2 *trail, int count)
{
    if (count < 2) return;
    for (int i = 1; i < count; i++)
    {
        float alpha = (float)i / count;  // 越旧越淡
        Color c = {0, 120, 220, (unsigned char)(alpha * 180)};
        DrawLineEx(trail[i - 1], trail[i], 2.0f, c);
    }
}

int main()
{
    InitWindow(WIDTH, HEIGHT, "Double Pendulum Simulation");
    SetTargetFPS(60);

    // 摆的初始状态
    Pendulum p = {PI / 2.0f, PI / 2.0f, 0.0f, 0.0f};
    Vector2 origin = {WIDTH / 2.0f, HEIGHT / 4.0f};

    // 轨迹缓冲区（环形）
    Vector2 trail[TRAIL_LENGTH] = {0};
    int trailHead = 0;
    int trailCount = 0;
    int trailSkip = 2;  // 每N帧记录一次轨迹

    bool paused = false;
    float timeScale = 1.0f;  // 时间倍率，上下方向键调整

    while (!WindowShouldClose())
    {
        float dt = GetFrameTime();
        if (dt > 1.0f / 30.0f) dt = 1.0f / 30.0f;  // 防止大帧跳跃

        // ---- 输入处理 ----
        if (IsKeyPressed(KEY_SPACE)) paused = !paused;
        if (IsKeyPressed(KEY_R))
        {
            // 重置
            p = (Pendulum){PI / 2.0f, PI / 2.0f, 0.0f, 0.0f};
            trailHead = 0;
            trailCount = 0;
            timeScale = 1.0f;
        }

        // 上下方向键调速
        if (IsKeyPressed(KEY_UP))   timeScale *= 1.5f;
        if (IsKeyPressed(KEY_DOWN)) timeScale /= 1.5f;
        if (timeScale < 0.1f)  timeScale = 0.1f;
        if (timeScale > 10.0f) timeScale = 10.0f;

        // ---- 物理更新 ----
        if (!paused)
        {
            float simDt = dt * timeScale;
            float subDt = simDt / SUBSTEPS;
            for (int s = 0; s < SUBSTEPS; s++)
                Step(&p, subDt);

            // 记录轨迹
            trailSkip--;
            if (trailSkip <= 0)
            {
                trailSkip = 2;
                Vector2 pos2 = GetBobPos(p.th2, GetBobPos(p.th1, origin, L1), L2);
                trail[trailHead] = pos2;
                trailHead = (trailHead + 1) % TRAIL_LENGTH;
                if (trailCount < TRAIL_LENGTH) trailCount++;
            }
        }

        // ---- 绘制 ----
        BeginDrawing();
        ClearBackground(RAYWHITE);

        // 标题
        DrawText("Double Pendulum Simulation", 10, 10, 20, DARKGRAY);

        // 信息栏
        DrawText(TextFormat("FPS: %d  |  Speed: %.1fx", GetFPS(), timeScale), 10, 35, 16, GRAY);
        DrawText("[UP/DOWN] Speed  [SPACE] Pause  [R] Reset", 10, 55, 16, GRAY);
        if (paused)
            DrawText("PAUSED", WIDTH / 2 - 40, HEIGHT / 2, 30, Fade(RED, 0.6f));

        // 轨迹
        if (trailCount > 0)
        {
            // 从 head 开始按时间顺序绘制（最旧在前）
            Vector2 ordered[TRAIL_LENGTH];
            int start = (trailCount < TRAIL_LENGTH) ? 0 : trailHead;
            for (int i = 0; i < trailCount; i++)
                ordered[i] = trail[(start + i) % TRAIL_LENGTH];
            DrawTrail(ordered, trailCount);
        }

        // 摆
        DrawPendulum(origin, p.th1, p.th2);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}