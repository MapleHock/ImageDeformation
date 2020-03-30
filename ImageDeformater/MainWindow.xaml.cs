using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Drawing;
using System.Text.RegularExpressions;
using System.IO;
using System.Diagnostics;

namespace ImageDeformater {
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : Window {
        private string currentOpenFilePath;
        private string TPScontrolFilePath;
        private class pointPairItem {
            public double controlX {get;set; }
            public double controlY {get;set; }
            public double targetX {get;set; }
            public double targetY {get;set; }
        }
        public MainWindow() {
            InitializeComponent();
            currentOpenFilePath = "";
        }

        private void OpenBtn_Click(object sender, RoutedEventArgs e) {
            Microsoft.Win32.OpenFileDialog fileDialog = new Microsoft.Win32.OpenFileDialog();
            fileDialog.DefaultExt = "*jpg";
            fileDialog.Filter = "image files |*.jpg;*.jpeg;*.png;*.bmp|All files (*.*)|*.*";
            string initPath = System.Reflection.Assembly.GetExecutingAssembly().Location.Replace("\\ImageDeformater.exe","");
            fileDialog.InitialDirectory = initPath;
            if (fileDialog.ShowDialog() == true) { // bool?转bool
                string fullfile = fileDialog.FileName;
                currentOpenFilePath = fullfile;
                Utilities.bitmapOriginal = new Bitmap(fullfile);
                Utilities.imgOriginal = Utilities.Bitmap2Matrix(Utilities.bitmapOriginal);
                showBitmap(Utilities.bitmapOriginal, OriginalImg);

                // recheck 变形方法 radiobutton，即运行一次自动参数生成
                if (RotationRdo.IsChecked == true) {
                    RotationRdo.IsChecked = false; RotationRdo.IsChecked = true;
                }
                if (DistortRdo.IsChecked == true) {
                    DistortRdo.IsChecked = false; DistortRdo.IsChecked = true;
                }
                if (TPSRdo.IsChecked == true) {
                    TPSRdo.IsChecked = false; TPSRdo.IsChecked = true;
                }
            }
        }

        private void showBitmap(Bitmap bitmap, System.Windows.Controls.Image ImageBox) {
            var bitmapImage = new BitmapImage();
            using (var memory = new System.IO.MemoryStream()) {
                bitmap.Save(memory, System.Drawing.Imaging.ImageFormat.Png);
                memory.Position = 0;
                bitmapImage.BeginInit();
                bitmapImage.StreamSource = memory;
                bitmapImage.CacheOption = BitmapCacheOption.OnLoad;
                bitmapImage.EndInit();
                bitmapImage.Freeze();
            }
            ImageBox.Source = bitmapImage;
        }

        private void RunBtn_Click(object sender, RoutedEventArgs e) {
            if (Utilities.imgOriginal == null) // 尚未加载图片
                return;

            InterpolationFunc interpolationFunc;
            switch (InterpolationCmb.SelectedIndex) {
                case 0:
                    interpolationFunc = Utilities.NearestNeighbour;
                    break;
                case 1:
                    interpolationFunc = Utilities.Bilinear;
                    break;
                case 2:
                    interpolationFunc = Utilities.BiCubic;
                    break;
                case 3:
                    interpolationFunc = Utilities.BiCubicFast;
                    break;
                default:
                    interpolationFunc = Utilities.BiCubicFast;
                    break;
            }
            Stopwatch calcStopWatch = new Stopwatch();
            // 根据radio button情况，调用不同变形函数，并且计算变形函数运算用时

            // 扭曲旋转：
            if (RotationRdo.IsChecked == true) {
                try {
                    double maxAngle = Convert.ToDouble(RotateMaxAngleText.Text) * Math.PI / 180; // 输入的角度转换为弧度制，然后进行后端函数调用
                    double radius = Convert.ToDouble(RotateRadiusText.Text);
                    double centerX = Convert.ToDouble(RotateCenterXText.Text);
                    double centerY = Convert.ToDouble(RotateCenterYText.Text);
                    if (IsParallelCbx.IsChecked == false) {
                        calcStopWatch.Start();
                        Utilities.imgChanged = Utilities.RotateTwistImg(Utilities.imgOriginal,
                                                interpolationFunc,
                                                maxAngle, radius,
                                                centerX, centerY);
                        calcStopWatch.Stop();
                    } else {
                        calcStopWatch.Start();
                        Utilities.imgChanged = Utilities.RotateTwistImgParallel(Utilities.imgOriginal,
                                                interpolationFunc,
                                                maxAngle, radius,
                                                centerX, centerY);
                        calcStopWatch.Stop();
                    }
                } catch (Exception) {
                    MessageBox.Show("填入了非法参数！");
                    return;
                }          
            }

            // 图像畸变：
            if (DistortRdo.IsChecked == true) {
                try {
                    double radius = Convert.ToDouble(DistortRadiusText.Text);
                    double centerX = Convert.ToDouble(DistortCenterXText.Text);
                    double centerY = Convert.ToDouble(DistortCenterYText.Text);
                    if (IsParallelCbx.IsChecked == false) {
                        calcStopWatch.Start();
                        Utilities.imgChanged = Utilities.DistortImg(Utilities.imgOriginal,
                                                interpolationFunc,
                                                radius,
                                                centerX, centerY,
                                                DistortIsPillowDistortCbx.IsChecked == true);
                        calcStopWatch.Stop();
                    } else {
                        calcStopWatch.Start();
                        Utilities.imgChanged = Utilities.DistortImgParallel(Utilities.imgOriginal,
                                                interpolationFunc,
                                                radius,
                                                centerX, centerY,
                                                DistortIsPillowDistortCbx.IsChecked == true);
                        calcStopWatch.Stop();
                    }
                    
                } catch (Exception) {
                    MessageBox.Show("填入了非法参数！");
                    return;
                }
                
            }

            // TPS
            if (TPSRdo.IsChecked == true) {
                int count = PointPairList.Items.Count;
                Utilities.TPSpointPair = new double[count, 4];
                count--;
                for (; count >= 0; count--) {
                    pointPairItem currentItem = (pointPairItem)PointPairList.Items[count];
                    Utilities.TPSpointPair[count, 0] = currentItem.controlX;
                    Utilities.TPSpointPair[count, 1] = currentItem.controlY;
                    Utilities.TPSpointPair[count, 2] = currentItem.targetX;
                    Utilities.TPSpointPair[count, 3] = currentItem.targetY;
                }
                if (IsParallelCbx.IsChecked == false) {
                    calcStopWatch.Start();
                    Utilities.imgChanged = Utilities.TPSImg(Utilities.imgOriginal, interpolationFunc, Utilities.TPSpointPair);
                    calcStopWatch.Stop();
                } else {
                    calcStopWatch.Start();
                    Utilities.imgChanged = Utilities.TPSImgParallel(Utilities.imgOriginal, interpolationFunc, Utilities.TPSpointPair);
                    calcStopWatch.Stop();
                }
                
            }

            pureRunTimeLbl.Content = "净运算用时" + calcStopWatch.ElapsedMilliseconds.ToString() + "ms";
            if (TrimRdo.IsChecked == true)
                Utilities.imgChanged = Utilities.TrimImg(Utilities.imgChanged);

            Utilities.bitmapChanged = Utilities.Matrix2Bitmap(Utilities.imgChanged);
            if (Utilities.bitmapChanged != null)
                showBitmap(Utilities.bitmapChanged, ChangedImg);
        }

        // 保存文件
        private void SaveBtn_Click(object sender, RoutedEventArgs e) {
            Microsoft.Win32.SaveFileDialog fileDialog = new Microsoft.Win32.SaveFileDialog();
            fileDialog.DefaultExt = "*jpg";
            fileDialog.Filter = "image files |*.jpg;*.jpeg;*.png;*.bmp|All files (*.*)|*.*";
            string[] currentDirPart = currentOpenFilePath.Split('\\');
            fileDialog.FileName = "Changed" + currentDirPart[currentDirPart.Length - 1];
            string initPath = System.Reflection.Assembly.GetExecutingAssembly().Location.Replace("\\ImageDeformater.exe", "");
            fileDialog.InitialDirectory = initPath;
            if (fileDialog.ShowDialog() == true) {
                string filePath = fileDialog.FileName;
                BitmapEncoder encoder = new PngBitmapEncoder();
                encoder.Frames.Add(BitmapFrame.Create((BitmapImage)ChangedImg.Source));
                using (var fileStream = new System.IO.FileStream(filePath, System.IO.FileMode.Create)) {
                    encoder.Save(fileStream);
                }
            }
        }

        // TPS控制点图像打开，若控制点图像有对应的存储关键点信息的.txt文件，一并读取到
        // 点对列表中
        private void ControlOpenBtn_Click(object sender, RoutedEventArgs e) {
            Microsoft.Win32.OpenFileDialog fileDialog = new Microsoft.Win32.OpenFileDialog();
            fileDialog.DefaultExt = "*jpg";
            fileDialog.Filter = "image files |*.jpg;*.jpeg;*.png;*.bmp|All files (*.*)|*.*";
            string initPath = System.Reflection.Assembly.GetExecutingAssembly().Location.Replace("\\ImageDeformater.exe", "");
            fileDialog.InitialDirectory = initPath;
            if (fileDialog.ShowDialog() != true)  // bool?转bool
                return;

            string fullfile = fileDialog.FileName;
            TPScontrolFilePath = fullfile;
            var bitmapControl = new Bitmap(fullfile);
            showBitmap(bitmapControl, TPSControlImg);

            // 重新生成点对列表
            PointPairList.Items.Clear();
            string controlInfo = Regex.Replace(TPScontrolFilePath, "\\..*", ".txt"); ;
            string targetInfo = Regex.Replace(currentOpenFilePath, "\\..*", ".txt");
            if (File.Exists(controlInfo) && File.Exists(targetInfo)) {
                StreamReader controlPointFile = new StreamReader(controlInfo, Encoding.Default);
                StreamReader targetPointFile = new StreamReader(targetInfo, Encoding.Default);
                string[] convertStrBuffer;
                double[] convertDoubleBuffer = new double[4];
                string controlBuffer = controlPointFile.ReadLine();
                string targetBuffer = targetPointFile.ReadLine();
                while (controlBuffer != null && targetBuffer != null) {
                    convertStrBuffer = controlBuffer.Split(" ".ToCharArray());
                    if (convertStrBuffer.Length < 2) // 防止不正确文件格式带来错误
                        break;

                    convertDoubleBuffer[0] = Convert.ToDouble(convertStrBuffer[1]);
                    convertDoubleBuffer[1] = Convert.ToDouble(convertStrBuffer[0]);

                    convertStrBuffer = targetBuffer.Split(" ".ToCharArray());
                    if (convertStrBuffer.Length < 2) // 防止不正确文件格式带来错误
                        break;

                    convertDoubleBuffer[2] = Convert.ToDouble(convertStrBuffer[1]);
                    convertDoubleBuffer[3] = Convert.ToDouble(convertStrBuffer[0]);
                    PointPairList.Items.Add(new pointPairItem
                                                { controlX = convertDoubleBuffer[0],
                                                  controlY = convertDoubleBuffer[1],
                                                  targetX = convertDoubleBuffer[2],
                                                  targetY = convertDoubleBuffer[3],
                                          });

                    controlBuffer = controlPointFile.ReadLine();
                    targetBuffer = targetPointFile.ReadLine();
                }
            }
        }

        private void MethodRdoGroup_Checked(object sender, RoutedEventArgs e) {
            var rdo = (RadioButton)sender;

            ChangedImg.Source = null; // 切换时不显示效果图
            TPSControlImg.Source = null;
            
            if (rdo == RotationRdo) {
                int width = 0;
                int height = 0;
                if (Utilities.bitmapOriginal != null) {
                    width = Utilities.bitmapOriginal.Width;
                    height = Utilities.bitmapOriginal.Height;
                }                   
                RotateMaxAngleText.IsEnabled = true;
                RotateMaxAngleText.Text = "30";
                RotateRadiusText.IsEnabled = true;
                RotateRadiusText.Text = Convert.ToString(Math.Min(width, height) / 4);
                RotateCenterXText.IsEnabled = true;
                RotateCenterXText.Text = Convert.ToString(height / 2);
                RotateCenterYText.IsEnabled = true;
                RotateCenterYText.Text = Convert.ToString(width / 2);

                DistortIsPillowDistortCbx.IsChecked = true;
                DistortIsPillowDistortCbx.IsEnabled = false;
                DistortRadiusText.Text = "";
                DistortRadiusText.IsEnabled = false;
                DistortCenterXText.Text = "";
                DistortCenterXText.IsEnabled = false;
                DistortCenterYText.Text = "";
                DistortCenterYText.IsEnabled = false;


                ControlOpenBtn.IsEnabled = false;
                PointPairList.Items.Clear();
                PointPairList.IsEnabled = false;
                TPSInsertNewBtn.IsEnabled = false;
                TPSNewPairText.Text = "";
                TPSNewPairText.IsEnabled = false;
                TPSControlImg.Source = null;
                TrimRdo.IsEnabled = false;
                return;
            }

            if(rdo == DistortRdo) {
                int width = 0;
                int height = 0;
                if (Utilities.bitmapOriginal != null) {
                    width = Utilities.bitmapOriginal.Width;
                    height = Utilities.bitmapOriginal.Height;
                }

                DistortIsPillowDistortCbx.IsEnabled = true;
                DistortIsPillowDistortCbx.IsChecked = true;
                DistortRadiusText.IsEnabled = true;
                DistortRadiusText.Text = Convert.ToString(Math.Min(width,height) / 2);
                DistortCenterXText.IsEnabled = true;
                DistortCenterXText.Text = Convert.ToString(height / 2);
                DistortCenterYText.IsEnabled = true;
                DistortCenterYText.Text = Convert.ToString(width / 2);


                RotateMaxAngleText.Text = "";
                RotateMaxAngleText.IsEnabled = false;
                RotateRadiusText.Text = "";
                RotateRadiusText.IsEnabled = false;
                RotateCenterXText.Text = "";
                RotateCenterXText.IsEnabled = false;
                RotateCenterYText.Text = "";
                RotateCenterYText.IsEnabled = false;

                
                ControlOpenBtn.IsEnabled = false;
                PointPairList.Items.Clear();
                PointPairList.IsEnabled = false;
                TPSInsertNewBtn.IsEnabled = false;
                TPSNewPairText.Text = "";
                TPSNewPairText.IsEnabled = false;
                TPSControlImg.Source = null;
                TrimRdo.IsEnabled = false;
                return;              
            }

            if (rdo == TPSRdo) {
                
                ControlOpenBtn.IsEnabled = true;
                PointPairList.IsEnabled = true;
                TPSInsertNewBtn.IsEnabled = true;
                TPSNewPairText.IsEnabled = true;
                TrimRdo.IsEnabled = true;

                RotateMaxAngleText.Text = "";
                RotateMaxAngleText.IsEnabled = false;
                RotateRadiusText.Text = "";
                RotateRadiusText.IsEnabled = false;
                RotateCenterXText.Text = "";
                RotateCenterXText.IsEnabled = false;
                RotateCenterYText.Text = "";
                RotateCenterYText.IsEnabled = false;

                DistortIsPillowDistortCbx.IsChecked = true;
                DistortIsPillowDistortCbx.IsEnabled = false;
                DistortRadiusText.Text = "";
                DistortRadiusText.IsEnabled = false;
                DistortCenterXText.Text = "";
                DistortCenterXText.IsEnabled = false;
                DistortCenterYText.Text = "";
                DistortCenterYText.IsEnabled = false;
                return;
            }
        }

        // 手动TPS接口，插入一组新的TPS点对
        private void TPSInsertNewBtn_Click(object sender, RoutedEventArgs e) {
            var strBuffer = TPSNewPairText.Text.Split(",".ToCharArray());
            if (strBuffer.Length == 4) {
                try {
                    PointPairList.Items.Add(new pointPairItem {
                        targetX = Convert.ToDouble(strBuffer[0]),
                        targetY = Convert.ToDouble(strBuffer[1]),
                        controlX = Convert.ToDouble(strBuffer[2]),
                        controlY = Convert.ToDouble(strBuffer[3]),
                    });
                } catch (Exception) {
                    MessageBox.Show("填入了非法参数！");
                    return;
                }
                
            }
            TPSNewPairText.Text = "";
        }

        // 辅助功能，点击图片时返回像素的坐标信息，用于配合手动TPS
        private void GetImgPixelInfo_MouseDown(object sender, MouseButtonEventArgs e) {
            var imageControl = (System.Windows.Controls.Image)sender;
            var bitmapImage = (BitmapImage)(imageControl.Source);
            int pixelMousePositionY = Convert.ToInt32(e.GetPosition(imageControl).X * bitmapImage.PixelWidth / imageControl.Width);
            int pixelMousePositionX = Convert.ToInt32(e.GetPosition(imageControl).Y * bitmapImage.PixelHeight / imageControl.Height);
            PosInfoTbx.Text = pixelMousePositionX.ToString() + "," + pixelMousePositionY.ToString();
        }
    }

    
}
