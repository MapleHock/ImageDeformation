using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Drawing;
using System.IO;

namespace ImageDeformater {
    static class Utilities {
        public static Bitmap bitmapOriginal;
        public static Bitmap bitmapChanged;
        public static byte[,,] imgOriginal;
        public static byte[,,] imgChanged;
        public static double[,] TPSpointPair;


        // 插值方法 参数格式byte[,,]img, double x, double y, int channel
        // 返回byte,为值类型返回。防止内存碎块,同时加快运行速度

        public static byte NearestNeighbour(byte[,,]img, double x, double y, int channel) {
            int nearestX = -1;
            int nearestY = -1;
            double roundX = Math.Round(x);
            double roundY = Math.Round(y);
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            // 超界点坐标插值设置,防止越界。
            if (roundX < double.Epsilon) 
                nearestX = 0;
            if (roundX >= height - double.Epsilon)
                nearestX = height - 1;
            if (roundY < double.Epsilon)
                nearestY = 0;
            if (roundY >= width - double.Epsilon)
                nearestY = 0;

            // 非超界点坐标插值设置
            if (nearestX == -1)
                nearestX = Convert.ToInt32(roundX);
            
            if (nearestY == -1)
                nearestY = Convert.ToInt32(roundY);

            return img[nearestX,nearestY,channel];
        }

        public static byte Bilinear(byte[,,]img, double x, double y, int channel) {
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            // 取整与判断边界值合并
            int top = Convert.ToInt32(Math.Min(Math.Max(Math.Floor(x),0),height - 1)); 
            int bot = Convert.ToInt32(Math.Min(Math.Max(Math.Ceiling(x), 0), height - 1));
            int left = Convert.ToInt32(Math.Min(Math.Max(Math.Floor(y), 0), height - 1));
            int right = Convert.ToInt32(Math.Min(Math.Max(Math.Ceiling(y), 0), height - 1));
            double portionTop;
            double portionLeft;
            if (top == bot) // 对应与超界情况，避免使用运算决定分比而导致误差
                portionTop = 1;
            else
                portionTop = bot - x;
            if (left == right)
                portionLeft = 1;
            else
                portionLeft = right - y;

            double grayValue = portionTop * (img[top, left, channel] * portionLeft + img[top, right, channel] * (1 - portionLeft))
                             + (1 - portionTop) * (img[bot, left, channel] * portionLeft + img[bot, right, channel] * (1 - portionLeft));
            
            return Convert.ToByte(grayValue);
        }

        public static byte BiCubic(byte[,,] img, double x, double y, int channel) {
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,] neighbourMatrix = new byte[4, 4];
            double[] cubicPortionX = new double[4];
            double[] cubicPortionY = new double[4];
            int baseX = Convert.ToInt32(Math.Min(Math.Max(Math.Floor(x), 0), height - 1));
            int baseY = Convert.ToInt32(Math.Min(Math.Max(Math.Floor(y), 0), width - 1));
            double u = x - baseX;
            double v = y - baseY;
            for (var i = 0; i < 4; i++) {
                int varX = Math.Min(Math.Max(baseX - 1 + i, 0), height - 1);
                for (var j = 0; j < 4; j++) {
                    int varY = Math.Min(Math.Max(baseY - 1 + j, 0), width - 1);
                    // 构造计算式ABC^T中的矩阵B
                    neighbourMatrix[i, j] = img[varX, varY, channel];

                    // 构造计算式ABC^T中的矩阵C
                    if (i == 0) {
                        cubicPortionY[j] = CubicBaseFunc(v + 1 - j);
                    }

                    // 构造计算式ABC^T中的矩阵A
                    if (j == 0) {
                        cubicPortionX[i] = CubicBaseFunc(u + 1 - i);
                    }
                }
            }
            double grayValue = 0;
            // 参照二次型或者双线性型的展开，最终结果为按行向量列表和列向量行标索引到系数矩阵，对应相乘再相加
            for (var i = 0; i < 4; i++) {
                for (var j = 0; j < 4; j++) {
                    grayValue += cubicPortionX[i] * cubicPortionY[j] * neighbourMatrix[i, j];
                }
            }
            grayValue = Math.Min(Math.Max(grayValue, 0), byte.MaxValue);
            return Convert.ToByte(grayValue);
        }

        private static double CubicBaseFunc(double x) {
            // 双三次插值权重函数
            double absX = Math.Abs(x);
            if (absX > 2)
                return 0;
            if (absX <= 1)               
                return 1 - absX * absX * (2 - absX); // 原始公式 (1 - 2 * absX * absX + absX * absX * absX)，按秦九韶算法修改后乘法少效率高

            return absX * (absX * (-absX + 5) - 8) + 4; // 原始公式 (4 - 8 * absX + 5 * absX * absX - absX * absX * absX)
        }

        public static byte BiCubicFast(byte[,,] img, double x, double y, int channel) {
            // 双三次插值快速版本，注意到代入的uv均为0,1之间的量，故而每一个代入权重计算函数的u, u + 1, u -1, u-2
            // 都是在权重函数固定区间，可以依次展开为关于u,v的多项式
            // 由此，直接在函数题内部，同uv的一二三次方直接内联计算权重向量即可
            // 这个修改思路减少了重复的次方运算，也减少了函数的调用，能提高插值速度
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,] neighbourMatrix = new byte[4, 4];
            double[] cubicPortionX = new double[4];
            double[] cubicPortionY = new double[4];
            int baseX = Convert.ToInt32(Math.Min(Math.Max(Math.Floor(x), 0), height - 1));
            int baseY = Convert.ToInt32(Math.Min(Math.Max(Math.Floor(y), 0), width - 1));
            double u = x - baseX;
            double v = y - baseY;
            double uSquare = u * u;
            double vSquare = v * v;
            double uCubic = u * uSquare;
            double vCubic = v * vSquare;
            // 对应于普通版本中的在循环中构造矩阵A、C
            cubicPortionX[0] = -u + 2 * uSquare - uCubic;
            cubicPortionX[1] = 1 - 2 * uSquare + uCubic;
            cubicPortionX[2] = u + uSquare - uCubic;
            cubicPortionX[3] = uCubic - uSquare;

            cubicPortionY[0] = -v + 2 * vSquare - vCubic;
            cubicPortionY[1] = 1 - 2 * vSquare + vCubic;
            cubicPortionY[2] = v + vSquare - vCubic;
            cubicPortionY[3] = vCubic - vSquare;

            double grayValue = 0;
            // 参照二次型或者双线性型的展开，最终结果为按行向量列表和列向量行标索引到系数矩阵，对应相乘再相加
            for (var i = 0; i < 4; i++) {
                for (var j = 0; j < 4; j++) {
                    grayValue += cubicPortionX[i] * cubicPortionY[j] 
                                        * img[Math.Min(Math.Max(baseX - 1 + i, 0), height - 1)
                                                        , Math.Min(Math.Max(baseY - 1 + j, 0), width - 1)
                                                        , channel];
                }
            }
            grayValue = Math.Min(Math.Max(grayValue, 0), byte.MaxValue);
            return Convert.ToByte(grayValue);
        }

        //------------------插值方法 End-----------------------------------------

        // 变形方法 参数格式byte[,,]img, InterpolationFunc interpolationFunc, 其他变形相关参数 返回值byte[,,]

        // maxAngle弧度制，在前端中，为了用户友好，规定写入角度值，并且由前端完成角度弧度转换这样的简单逻辑
        public static byte[,,] RotateTwistImg(byte[,,] img, InterpolationFunc interpolationFunc, 
                                double maxAngle, double radius, 
                                double centerX, double centerY) {
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,,] imgTwisted = new byte[height, width, img.GetLength(2)];  
            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    // 按照变化公式，计算目标图(i,j)应到原图哪个位置取像素信息
                    double distance = Math.Sqrt((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY));
                    if (distance > radius) {
                        for (int channel = 0; channel < 3; channel++) {
                            imgTwisted[i, j, channel] = img[i, j, channel];
                        }
                        continue;
                    }
                    // 原始公式 angle = maxAngle * (radius - distance) / radius， 化简防止radius distance差距过大大数吃小数
                    double angle = maxAngle * (1 - distance / radius) ; 
                    double x = Math.Cos(angle) * (i - centerX) - Math.Sin(angle) * (j - centerY) + centerX;
                    double y = Math.Sin(angle) * (i - centerX) + Math.Cos(angle) * (j - centerY) + centerY;
                    for (int channel = 0; channel < 3; channel++) {
                        imgTwisted[i, j, channel] = interpolationFunc(img, x, y, channel);
                    }
                }
            }
            return imgTwisted;
        }


        public static byte[,,] DistortImg(byte[,,] img, InterpolationFunc interpolationFunc,
                               double radius,
                               double centerX, double centerY,
                               bool isAdjustPillowDistort) {
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,,] imgDistorted = new byte[height, width, img.GetLength(2)];
            if (isAdjustPillowDistort) { // 校正枕形畸变，即对图像做桶形畸变
                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++) {
                        double distance = Math.Sqrt((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY));
                        if (distance > radius) {
                            continue; // 投影球外的点，置为黑 0,0,0，即默认值
                        }
                        double k;
                        if (distance / radius < 0.001)
                            k = 1; // k 极小时
                        else
                            k = radius / distance * Math.Asin(distance / radius);

                        // 为了精确快速调整计算顺序，把原始公式在注释给出
                        double x = k * i + (1 - k) * centerX; // 原始公式 k * (i - centerX) + centerX
                        double y = k * j + (1 - k) * centerY; // 原始公式 k * (i - centerX) + centerX
                        if (x >= height || x < 0 || y >= width || y < 0)
                            continue;
                        for (int channel = 0; channel < 3; channel++) {
                            imgDistorted[i, j, channel] = interpolationFunc(img, x, y, channel);
                        }
                    }
                }
            } else { // 校正桶形畸变，即对图像做枕形畸变
                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++) {
                        double distance = Math.Sqrt((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY));
                        if (distance > Math.PI / 2 *  radius) {
                            continue; // 投影球外的点，置为黑 0,0,0，即默认值
                        }
                        double k;
                        if (distance / radius < 0.001)
                            k = 1; // k 极小时
                        else
                            k = Math.Sin(distance / radius) * radius / distance;

                        // 为了精确快速调整计算顺序，把原始公式在注释给出
                        double x = k * i + (1 - k) * centerX; // 原始公式 k * (i - centerX) + centerX
                        double y = k * j + (1 - k) * centerY; // 原始公式 k * (i - centerX) + centerX
                        if (x >= height || x < 0 || y >= width || y < 0)
                            continue;
                        for (int channel = 0; channel < 3; channel++) {
                            imgDistorted[i, j, channel] = interpolationFunc(img, x, y, channel);
                        }
                    }
                }
            }
            return imgDistorted;
        }

        // pointPair[i,0] ~ [i,3] 分别为 控制点x1 控制点y1 目标点x2 目标点y2，控制点：想变成的参考样貌；目标点：待变形图上的坐标
        public static byte[,,] TPSImg(byte[,,] img, InterpolationFunc interpolationFunc,
                                double[,] pointPair) {
            if (pointPair.GetLength(1) != 4)
                return null; // Error pointPair

            int pointNumber = pointPair.GetLength(0);
            // 构造矩阵Y
            double[,] Y = new double[pointNumber + 3, 2];
            for (var i = 0; i < pointNumber; i++) {
                Y[i, 0] = pointPair[i, 2];
                Y[i, 1] = pointPair[i, 3];
            }

            // 构造矩阵L，注意到L为对称阵，赋值时，循环体可以在第二重循环时减半循环量
            double[,] L = new double[pointNumber + 3, pointNumber + 3];
            for (var j = 0; j < pointNumber; j++) {
                L[pointNumber, j] = 1;
                L[j, pointNumber] = 1;
            }
            for (var i = pointNumber + 1; i < pointNumber + 3; i++) {
                for (var j = 0; j < pointNumber; j++) {
                    L[i, j] = pointPair[j, i - pointNumber - 1];
                    L[j, i] = L[i, j];
                }
            }
            for (var i = 0; i < pointNumber; i++) {
                for (var j = 0; j <= i; j++) {
                    if (j == i) {
                        L[i, j] = 0;
                        break;
                    }
                    L[i, j] = RadialbasisByPoint(pointPair[i, 0], pointPair[i, 1], pointPair[j, 0], pointPair[j, 1]);
                    L[j, i] = L[i, j];
                }
            }
            // 传入高斯消元函数的为引用，会被函数体修改矩阵值，原本应使用
            // double[,] coe = GaussianElimination(L.Clone(), Y.Clone());
            // 但LY不再使用，故为节约耗时，不clone
            double[,] coe = GaussianElimination(L, Y);
            if (coe == null)
                return null; // 方程无解，返回空图

            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,,] imgTPSed = new byte[height, width, 3];
            double x;
            double y;
            for (var i = 0; i < height; i++) {
                for (var j = 0; j < width; j++) {
                    // 先计算 wU之外的项目
                    x = coe[pointNumber, 0] + coe[pointNumber + 1, 0] * i + coe[pointNumber + 2, 0] * j;
                    y = coe[pointNumber, 1] + coe[pointNumber + 1, 1] * i + coe[pointNumber + 2, 1] * j;
                    // 分别计算各wU的结果，按正负号存入数组中，累加时，
                    // 保证xy和新作为加项的wU正负异号
                    // 则整体加法流程充分地让正负相消累计，每一步得到的数的绝对值较小，
                    // 防止double型在位数提高时新引入舍入误差
                    double[,] weightedUArray = new double[pointNumber, 4];
                    double u;
                    double weightedUX;
                    double weightedUY;
                    int posUXcount = 0;
                    int negUXcount = 0;
                    int posUYcount = 0;
                    int negUYcount = 0;
                    // 计算wU并按正负号存储
                    for (var wCount = 0; wCount < pointNumber; wCount++) {
                        u = RadialbasisByPoint(i, j, pointPair[wCount, 0], pointPair[wCount, 1]);
                        weightedUX = u * coe[wCount, 0];
                        weightedUY = u * coe[wCount, 1];
                        if (weightedUX > 0) {
                            weightedUArray[posUXcount, 0] = weightedUX;
                            posUXcount++;
                        } else {
                            weightedUArray[negUXcount, 1] = weightedUX;
                            negUXcount++;
                        }

                        if (weightedUY > 0) {
                            weightedUArray[posUYcount, 2] = weightedUY;
                            posUYcount++;
                        } else {
                            weightedUArray[negUYcount, 3] = weightedUY;
                            negUYcount++;
                        }
                    }
                    posUXcount--; posUYcount--; negUXcount--; negUYcount--;
                    // 分正负，尽可能异号累加得到x
                    while (posUXcount >= 0 && negUXcount >= 0) {
                        if (x > 0) {
                            x += weightedUArray[negUXcount, 1];
                            negUXcount--;
                        } else {
                            x += weightedUArray[posUXcount, 0];
                            posUXcount--;
                        }
                    }
                    while (posUXcount >= 0) {
                        x += weightedUArray[posUXcount, 0];
                        posUXcount--;
                    }
                    while (negUXcount >= 0) {
                        x += weightedUArray[negUXcount, 1];
                        negUXcount--;
                    }

                    // 分正负，尽可能异号累加得到y
                    while (posUYcount >= 0 && negUYcount >= 0) {
                        if (y > 0) {
                            y += weightedUArray[negUYcount, 3];
                            negUYcount--;
                        } else {
                            y += weightedUArray[posUYcount, 2];
                            posUYcount--;
                        }
                    }
                    while (posUYcount >= 0) {
                        y += weightedUArray[posUYcount, 2];
                        posUYcount--;
                    }
                    while (negUYcount >= 0) {
                        y += weightedUArray[negUYcount, 3];
                        negUYcount--;
                    }

                    // 得到对应坐标点，到原图中查找像素信息
                    if (x < 0 || x >= height || y < 0 || y > width)
                        continue;
                    for (var channel = 0; channel < 3; channel++) {
                        imgTPSed[i, j, channel] = interpolationFunc(img, x, y, channel);
                    }
                }
            }
            return imgTPSed;
        }

        // 辅助函数，解线性方程组，输入点组的径向基，会改变矩阵值，若不想被引用修改可通过传入clone解决
        // 返回null无解 
        private static double[,] GaussianElimination(double[,] A, double[,] b) {
            if (A.GetLength(0) != A.GetLength(1) || A.GetLength(0) != b.GetLength(0))
                return null; // error rank

            int unknownNum = b.GetLength(0);
            int groupNum = b.GetLength(1);
            for (var i = 0; i < unknownNum; i++) {
                if (Math.Abs(A[i,i]) < Math.Pow(10, -6)) {
                    var iBelow = i + 1;
                    for (; iBelow < unknownNum; iBelow++) {
                        // 列主元较小，交换行
                        if (Math.Abs(A[iBelow,i]) > Math.Pow(10, -6)) {
                            double temp;
                            for (var j = i; j < unknownNum; j++) {
                                temp = A[iBelow, j];
                                A[iBelow, j] = A[i,j];
                                A[i, j] = temp;
                            }
                            for (var j = 0; j < groupNum; j++) {
                                temp = b[iBelow, j];
                                b[iBelow, j] = b[i, j];
                                b[i, j] = temp;
                            }
                            break;
                        }
                    }
                    if (iBelow == unknownNum)
                        return null;
                }

                for (var iBelow = i + 1; iBelow < unknownNum; iBelow++) {
                    double eCoe = - A[iBelow, i] / A[i, i];
                    A[iBelow, i] = 0;
                    for (var j = i + 1; j < unknownNum; j++) {
                        A[iBelow, j] += eCoe * A[i, j];
                    }
                    for (var j = 0; j < groupNum; j++) {
                        b[iBelow, j] += eCoe * b[i, j];
                    }
                }
            }

            for (var i = unknownNum - 1; i >= 0; i--) {
                for (var j = 0; j < groupNum; j++) {
                    b[i, j] /= A[i, i];
                    // 数学上还应有
                    // A[i, i] = 1;
                    // 但为了节约计算用时，不做这步
                }
                for (var iAbove = i - 1; iAbove >=0; iAbove--) {                  
                    for (var j = 0; j < groupNum; j++) {
                        b[iAbove, j] -= A[iAbove, i] * b[i, j];
                    }
                    // 数学上还应有
                    // A[iAbove, i] = 0;
                    // 但为了节约计算用时，不做这步
                }
            }

            return b;
        }

        private static double RadialbasisByPoint(double x1, double y1, double x2, double y2) {
            // 径向基函数
            double distanceSquare = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
            if (distanceSquare < Math.Pow(10, -6))
                return 0;
            return distanceSquare * Math.Log(distanceSquare);
        }

        //------------------变形方法和其辅助函数 End-----------------------------------------

        // 三种变形的并行版本
        // 变换思路和代码基本同前，只将为目标图每个格点赋像素值的步骤进行并行，即每个格点同时
        // 计算其在原图中对应的坐标，并查找像素信息
        // 故只将串行的两重for，转换为并行的循环

        public static byte[,,] RotateTwistImgParallel(byte[,,] img, InterpolationFunc interpolationFunc,
                                double maxAngle, double radius,
                                double centerX, double centerY) {
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,,] imgTwisted = new byte[height, width, img.GetLength(2)];
            System.Threading.Tasks.Parallel.For(0, height * width, index => {
                int i = index / width;
                int j = index % width;
                double distance = Math.Sqrt((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY));
                if (distance > radius) {
                    for (int channel = 0; channel < 3; channel++) {
                        imgTwisted[i, j, channel] = img[i, j, channel];
                    }
                    return;
                }
                double angle = maxAngle * (radius - distance) / radius;
                double x = Math.Cos(angle) * (i - centerX) - Math.Sin(angle) * (j - centerY) + centerX;
                double y = Math.Sin(angle) * (i - centerX) + Math.Cos(angle) * (j - centerY) + centerY;
                for (int channel = 0; channel < 3; channel++) {
                    imgTwisted[i, j, channel] = interpolationFunc(img, x, y, channel);
                }
            });
            return imgTwisted;
        }

        public static byte[,,] DistortImgParallel(byte[,,] img, InterpolationFunc interpolationFunc,
                               double radius,
                               double centerX, double centerY,
                               bool isAdjustPillowDistort) {
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,,] imgDistorted = new byte[height, width, img.GetLength(2)];
            if (isAdjustPillowDistort) { // 校正枕形畸变，即对图像做桶形畸变
                System.Threading.Tasks.Parallel.For(0, height * width, index => {
                    int i = index / width;
                    int j = index % width;
                    double distance = Math.Sqrt((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY));
                    if (distance > radius) {
                        return; // 投影球外的点，置为黑 0,0,0，即默认值
                    }
                    double k;
                    if (distance / radius < 0.001)
                        k = 1; // k 极小时
                    else
                        k = radius / distance * Math.Asin(distance / radius);

                    // 为了精确快速调整计算顺序，把原始公式在注释给出
                    double x = k * i + (1 - k) * centerX; // 原始公式 k * (i - centerX) + centerX
                    double y = k * j + (1 - k) * centerY; // 原始公式 k * (i - centerX) + centerX
                    if (x >= height || x < 0 || y >= width || y < 0)
                        return;
                    for (int channel = 0; channel < 3; channel++) {
                        imgDistorted[i, j, channel] = interpolationFunc(img, x, y, channel);
                    }
                });
                
            } else { // 校正桶形畸变，即对图像做枕形畸变
                System.Threading.Tasks.Parallel.For(0, height * width, index => {
                    int i = index / width;
                    int j = index % width;
                    double distance = Math.Sqrt((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY));
                    if (distance > Math.PI / 2 * radius) {
                        return; // 投影球外的点，置为黑 0,0,0，即默认值
                    }
                    double k;
                    if (distance / radius < 0.001)
                        k = 1; // k 极小时
                    else
                        k = Math.Sin(distance / radius) * radius / distance;

                    // 为了精确快速调整计算顺序，把原始公式在注释给出
                    double x = k * i + (1 - k) * centerX; // 原始公式 k * (i - centerX) + centerX
                    double y = k * j + (1 - k) * centerY; // 原始公式 k * (i - centerX) + centerX
                    if (x >= height || x < 0 || y >= width || y < 0)
                        return;
                    for (int channel = 0; channel < 3; channel++) {
                        imgDistorted[i, j, channel] = interpolationFunc(img, x, y, channel);
                    }

                });
            }
            return imgDistorted;
        }

        public static byte[,,] TPSImgParallel(byte[,,] img, InterpolationFunc interpolationFunc,
                                double[,] pointPair) {
            // 置矩阵，高斯消元解系数同串行版本
            if (pointPair.GetLength(1) != 4)
                return null; // Error pointPair

            int pointNumber = pointPair.GetLength(0);
            double[,] Y = new double[pointNumber + 3, 2];
            for (var i = 0; i < pointNumber; i++) {
                Y[i, 0] = pointPair[i, 2];
                Y[i, 1] = pointPair[i, 3];
            }

            double[,] L = new double[pointNumber + 3, pointNumber + 3];
            for (var j = 0; j < pointNumber; j++) {
                L[pointNumber, j] = 1;
                L[j, pointNumber] = 1;
            }
            for (var i = pointNumber + 1; i < pointNumber + 3; i++) {
                for (var j = 0; j < pointNumber; j++) {
                    L[i, j] = pointPair[j, i - pointNumber - 1];
                    L[j, i] = L[i, j];
                }
            }
            for (var i = 0; i < pointNumber; i++) {
                for (var j = 0; j <= i; j++) {
                    if (j == i) {
                        L[i, j] = 0;
                        break;
                    }
                    L[i, j] = RadialbasisByPoint(pointPair[i, 0], pointPair[i, 1], pointPair[j, 0], pointPair[j, 1]);
                    L[j, i] = L[i, j];
                }
            }
            // 传入高斯消元函数的为引用，会被函数体修改矩阵值，原本应使用
            // double[,] coe = GaussianElimination(L.Clone(), Y.Clone());
            // 但LY不再使用，故为节约耗时，不clone
            double[,] coe = GaussianElimination(L, Y);
            if (coe == null)
                return null;

            // 以下开始并行
            int height = img.GetLength(0);
            int width = img.GetLength(1);
            byte[,,] imgTPSed = new byte[height, width, 3];

            System.Threading.Tasks.Parallel.For(0, height * width, index => {
                double x;
                double y;
                int i = index / width;
                int j = index % width;
                x = coe[pointNumber, 0] + coe[pointNumber + 1, 0] * i + coe[pointNumber + 2, 0] * j;
                y = coe[pointNumber, 1] + coe[pointNumber + 1, 1] * i + coe[pointNumber + 2, 1] * j;
                double[,] weightedUArray = new double[pointNumber, 4]; // 分正负相消累计，防止大数吃小数
                double u;
                double weightedUX;
                double weightedUY;
                int posUXcount = 0;
                int negUXcount = 0;
                int posUYcount = 0;
                int negUYcount = 0;
                for (var wCount = 0; wCount < pointNumber; wCount++) {
                    u = RadialbasisByPoint(i, j, pointPair[wCount, 0], pointPair[wCount, 1]);
                    weightedUX = u * coe[wCount, 0];
                    weightedUY = u * coe[wCount, 1];
                    if (weightedUX > 0) {
                        weightedUArray[posUXcount, 0] = weightedUX;
                        posUXcount++;
                    } else {
                        weightedUArray[negUXcount, 1] = weightedUX;
                        negUXcount++;
                    }

                    if (weightedUY > 0) {
                        weightedUArray[posUYcount, 2] = weightedUY;
                        posUYcount++;
                    } else {
                        weightedUArray[negUYcount, 3] = weightedUY;
                        negUYcount++;
                    }
                }
                posUXcount--; posUYcount--; negUXcount--; negUYcount--;
                while (posUXcount >= 0 && negUXcount >= 0) {
                    if (x > 0) {
                        x += weightedUArray[negUXcount, 1];
                        negUXcount--;
                    } else {
                        x += weightedUArray[posUXcount, 0];
                        posUXcount--;
                    }
                }
                while (posUXcount >= 0) {
                    x += weightedUArray[posUXcount, 0];
                    posUXcount--;
                }
                while (negUXcount >= 0) {
                    x += weightedUArray[negUXcount, 1];
                    negUXcount--;
                }

                while (posUYcount >= 0 && negUYcount >= 0) {
                    if (y > 0) {
                        y += weightedUArray[negUYcount, 3];
                        negUYcount--;
                    } else {
                        y += weightedUArray[posUYcount, 2];
                        posUYcount--;
                    }
                }
                while (posUYcount >= 0) {
                    y += weightedUArray[posUYcount, 2];
                    posUYcount--;
                }
                while (negUYcount >= 0) {
                    y += weightedUArray[negUYcount, 3];
                    negUYcount--;
                }

                if (x < 0 || x >= height || y < 0 || y > width)
                    return;
                for (var channel = 0; channel < 3; channel++) {
                    imgTPSed[i, j, channel] = interpolationFunc(img, x, y, channel);
                }
            });


            return imgTPSed;
        }

        // 转换函数, 从bitmap坐标系x正向向右，y正向向下到matlab坐标系x正向下，y正向右（x行指标，y列指标）
        public static byte[,,] Bitmap2Matrix(Bitmap bitmap) {
            
            int height = bitmap.Height;
            int width = bitmap.Width;
            byte[,,] matrix = new byte[height, width, 3];
            for (var i = 0; i < height; i++) {
                for (var j = 0; j < width; j++) {
                    matrix[i, j, 0] = bitmap.GetPixel(j, i).R;
                    matrix[i, j, 1] = bitmap.GetPixel(j, i).G;
                    matrix[i, j, 2] = bitmap.GetPixel(j, i).B;
                }
            }

            return matrix;

        }

        public static Bitmap Matrix2Bitmap(byte[,,] matrix) {
            if (matrix == null) 
                return null;
           
            int height = matrix.GetLength(0);
            int width = matrix.GetLength(1);
            var bitmap = new Bitmap(width, height);
            for (var i = 0; i < height; i++) {
                for (var j = 0; j < width; j++) {
                    var temp = System.Drawing.Color.FromArgb(
                                    matrix[i, j, 0],
                                    matrix[i, j, 1],
                                    matrix[i, j, 2]);
                    bitmap.SetPixel(j, i, System.Drawing.Color.FromArgb(
                                    matrix[i, j, 0],
                                    matrix[i, j, 1],
                                    matrix[i, j, 2]));
                }
            }
            return bitmap;
        }

        // 辅助函数，用于裁剪图像黑边
        public static byte[,,] TrimImg(byte[,,] img) {
            int originalHeight = img.GetLength(0);
            int originalWidth = img.GetLength(1);

            int startTop = 0;
            for (;startTop < originalHeight;startTop++) {
                int colCount;
                for (colCount = 0; colCount < originalWidth; colCount++) {
                    if (img[startTop, colCount, 0] != 0 || img[startTop, colCount, 1] != 0 || img[startTop, colCount, 2] != 0)
                        break;
                }
                if (colCount != originalWidth)
                    break;
            }

            int startLeft = 0;
            for (;startLeft < originalWidth; startLeft++) {
                int rowCount;
                for (rowCount = 0; rowCount < originalHeight; rowCount++) {
                    if (img[rowCount, startLeft, 0] != 0 || img[rowCount, startLeft, 1] != 0 || img[rowCount, startLeft, 2] != 0)
                        break;
                }
                if (rowCount != originalHeight)
                    break;
            }

            int startBot = originalHeight - 1;
            for (;startBot >= 0; startBot--) {
                int colCount;
                for (colCount = 0; colCount < originalWidth; colCount++) {
                    if (img[startBot, colCount, 0] != 0 || img[startBot, colCount, 1] != 0 || img[startBot, colCount, 2] != 0)
                        break;
                }
                if (colCount != originalWidth)
                    break;
            }

            int startRight = originalWidth - 1;
            for (;startRight >= 0; startRight--) {
                int rowCount;
                for (rowCount = 0; rowCount < originalHeight; rowCount++) {
                    if (img[rowCount, startRight, 0] != 0 || img[rowCount, startRight, 1] != 0 || img[rowCount, startRight, 2] != 0)
                        break;
                }
                if (rowCount != originalHeight)
                    break;
            }

            int rows = startBot - startTop + 1;
            int cols = startRight - startLeft + 1;

            if (rows < 0 || cols < 0)
                return null;

            byte[,,] imgTrimed = new byte[rows, cols, 3];

            System.Threading.Tasks.Parallel.For(0, rows * cols, index => {
                int i = index / cols;
                int j = index % cols;
                for (int channel = 0; channel < 3; channel++) {
                    imgTrimed[i, j, channel] = img[i + startTop, j + startLeft, channel];
                }
            });

            return imgTrimed;
        }
    }
    // 插值函数委托，用作变形函数参数，让变形和插值两块解耦
    public delegate byte InterpolationFunc(byte[,,] img, double x, double y, int channel);
}
