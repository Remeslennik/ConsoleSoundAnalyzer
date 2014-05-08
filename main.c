#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "portaudio.h"
//файлы собранной библиотеки скачаны в конце статьи http://howtodoit.com.ua/kak-vosproizvesti-nepreryivnyie-zvukovyie-dannyie/

#define SAMPLE_RATE (44100)
#define FRAMES_PER_BUFFER (4096)
#define POW_2 (12) // для фурье берём 2 в степени 12 = 4096 семплов в буфере
#define WINDOWMA (22) // скользящая средняя по 22 семплам

/* Select sample format. */
#define PA_SAMPLE_TYPE paFloat32
#define SAMPLE_SIZE (4)// по 4 байта на семпл, 1=максимальное значение, -1=минимальное
#define SAMPLE_SILENCE (0.0f) // 0 - тишина
#define CLEAR(a) memset( (a), 0, FRAMES_PER_BUFFER * NUM_CHANNELS * SAMPLE_SIZE ) // обнуление массива
#define PRINTF_S_FORMAT "%.8f"

#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...
#define  FT_DIRECT        -1    // Direct transform. Прямое преобразование фурье
#define  FT_INVERSE        1    // Inverse transform.


int FFT(float *Rdat, float *Idat, int N, int LogN, int Ft_Flag); // преобразование фурье, код отсюда http://ru.wikipedia.org/wiki/Быстрое_преобразование_Фурье
void WIN_FUN(float *Index, int n); // коэффициенты оконной функции, реализация разных окон http://dmilvdv.narod.ru/SpeechSynthesis/window_function.html
void SimpleMA(float *buf_1, int buf_size, int window); // скользящее среднее, среднее арифметическое из WINDOWMA значений
void VerPar(float x1,float y1,float x2,float y2,float x3,float y3); // поиск вершины параболы по трем точкам
float GetFreq(float *Rdat, float *Idat, int N); // находим основную частоту. самая проблемная функция
float X0, Y0; // координаты вершины параболы

int main()
{
    PaStreamParameters inputParameters;
    PaStream *stream = NULL;

    float *sampleBlock;
    int i;
    int numBytes,NoteMIDI,delta;
    float RealNote,RealFreq;
    char notes[12][5]={"do","do#","re","re#","mi","fa","fa#","sol","sol#","la","la#","si"};
    float Im[FRAMES_PER_BUFFER]; // мнимая часть
    float WinFun[FRAMES_PER_BUFFER]; // оконная функция
    WIN_FUN(WinFun,FRAMES_PER_BUFFER); // вычисляем коэффициенты оконой функции

    numBytes = FRAMES_PER_BUFFER * NUM_CHANNELS * SAMPLE_SIZE ;
    sampleBlock = (float *) malloc( numBytes );
    CLEAR(sampleBlock);

    // инициализация портаудио
    Pa_Initialize();
    inputParameters.device = Pa_GetDefaultInputDevice();
    inputParameters.channelCount = NUM_CHANNELS;
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency ;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    // открытие потока
    Pa_OpenStream(
        &stream,
        &inputParameters,
        0,
        SAMPLE_RATE,
        FRAMES_PER_BUFFER,
        paClipOff,
        NULL,
        NULL );

    // старт потока
    Pa_StartStream( stream );
    printf("**************************************************\n");
    printf("               READING STREAM.\n");
    printf("**************************************************\n");
    
    while(1) // основой цикл, без выхода
    {
        Pa_ReadStream( stream, sampleBlock, FRAMES_PER_BUFFER ); // захват звука
        CLEAR( Im ); // обнуляем мнимую часть
        SimpleMA( sampleBlock, FRAMES_PER_BUFFER, WINDOWMA); // фильтр скользящее среднее
        for (i=0; i<FRAMES_PER_BUFFER; ++i) sampleBlock[i]*=WinFun[i];// применяем оконную функцию
        FFT(sampleBlock, Im, FRAMES_PER_BUFFER, POW_2, -1); // отправляем на фурье
        RealFreq=GetFreq(sampleBlock, Im, FRAMES_PER_BUFFER); // ищем частоту
        
        if (Y0 > 50 ) // если в спектре есть сингал
        {
            RealNote = 21 + (log( RealFreq / 27.5 )) / (log( pow (2, 1/12.))); // преобразуем в дробную midi ноту
            RealNote = floor(RealNote+0.5); // округляем
            delta=floor ((RealNote -NoteMIDI)*100)+0.5; // отклонение от ближайщей идеальной ноты, в музыкальных центах

            //выводим в консоль
            printf("Pow=%.1f ",Y0); // мощность максимальной частоты
            printf("Freq=%.3f ",RealFreq); // основная частота
            printf("note=%i ",NoteMIDI); // таблица соответствия - https://raw.github.com/xintrea/mytetra_syncro/master/base/1343982148ot0wjdzyil/image17192.png
            printf("%s ",notes[NoteMIDI%12]);  // название ноты
            printf("delta=%i ",delta); // отклонение в центах
            printf("\n");
            
        }
        
        
    CLEAR( sampleBlock );
    free( sampleBlock );
    Pa_StopStream( stream );
    Pa_Terminate();
    return 0;
}

void VerPar(float x1,float y1,float x2,float y2,float x3,float y3)
{
    // определяем координаты вершины параболы, выводим в глобал переменные
    float A,B,C;
    A=(y3-( (x3*(y2-y1)+x2*y1-x1*y2) /  (x2-x1)))   /(x3*(x3-x1-x2)+x1*x2);
    B= (y2-y1)/(x2-x1) - A*(x1+x2);
    C= ( (x2*y1-x1*y2) / (x2-x1) )+A*x1*x2;
    X0= -(B/(2*A));
    Y0= - ( (B*B-4*A*C)  /   (4*A) );
    return;
}

void SimpleMA(float *buf_1, int buf_size, int window)
{
    int i;
    float buf_2[buf_size];
    float sum=0;
    // начальные несколько семплов
    for (i=0; i<window; i++)
    {
        sum+=buf_1[i];
        buf_2[i]=sum/(i+1);
    }
    
    //основной буфер
    for (i=window; i<buf_size; i++)
    {
        sum-=buf_1[i-window]; // рекурентно, вычитаем старое значение
        sum+=buf_1[i]; // прибавляем новое
        buf_2[i]=sum/window;
    }

    for (i=0; i<buf_size; i++) buf_1[i]=buf_2[i];
    return;
}

void WIN_FUN(float *Index, int n)
{
    // формулы взяты с сайта http://dmilvdv.narod.ru/SpeechSynthesis/index.html?window_function.html
    // вычисляем коэффициенты оконной функции
    int i;
    double a;
    for (i=0; i<n ; i++)
    {
         a = M_PI / n * (1.0 + 2.0 * i);
         Index[i]=0.5 * (1.0 - 0.16 - cos(a) + 0.16 * cos(2.0 * a));
    }
    return;
}


float GetFreq (float *Rdat, float *Idat, int N)
{
    float max_ampl;
    int i, max_i;
    int i_start=60*FRAMES_PER_BUFFER/SAMPLE_RATE;// нижняя частота, примерный диапазон гитары
    int i_last=700*FRAMES_PER_BUFFER/SAMPLE_RATE;// верхняя частота
    float PowFreq[FRAMES_PER_BUFFER];
    float Freq_1=0, Freq_2=0, RealFreq=0;
    int i_2;
    float ampl_2;
    // определяем частоту с максимальным вектором
    max_ampl=0;
    max_i=0;
    for(i=i_start; i<i_last; i++)
    {
        PowFreq[i]=Rdat[i]*Rdat[i]+Idat[i]*Idat[i];
        if (PowFreq[i]>max_ampl)
        {
            max_ampl=PowFreq[i]; // максимальное значение
            max_i=i;             // и его индекс
        }
    }
    // интерполяция, определяем вершину параболы по 3 точкам для точного значения частоты
    VerPar((float)max_i-1,
           sqrt(  PowFreq[max_i-1] ),
           (float)max_i,
           sqrt(  max_ampl ),
           (float)max_i+1,
           sqrt(  PowFreq[max_i+1] ));

    Freq_1=X0*SAMPLE_RATE/FRAMES_PER_BUFFER; // вершина параболы в частоту
    Freq_2=Freq_1 / 2; // частота на октаву ниже для определения басов

    // ищем индексы и мощности гармонических частот
    i_2=floor((X0/2)+0.5); // ближайшая частота к басовой
    ampl_2=sqrt( PowFreq[i_2] ); // её значение
    
    // анализиурем по гармоникам
    RealFreq=Freq_1;
    if(Freq_1<1100)
    {
        // если основная частота тише второй гармоники
        if (ampl_2 > 30  && Freq_2>80)
        {
            RealFreq=Freq_2;
        }
    }
    return (RealFreq);
}


int  FFT(float *Rdat, float *Idat, int N, int LogN, int Ft_Flag)
{
    // parameters error check:
    if((Rdat == NULL) || (Idat == NULL))                  return 0;
    if((N > 16384) || (N < 1))                            return 0;
    if(!NUMBER_IS_2_POW_K(N))                             return 0;
    if((LogN < 2) || (LogN > 14))                         return 0;
    if((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE)) return 0;

    register int  i, j, n, k, io, ie, in, nn;
    float         ru, iu, rtp, itp, rtq, itq, rw, iw, sr;

    static const float Rcoef[14] =
    {
        -1.0000000000000000F,  0.0000000000000000F,  0.7071067811865475F,
        0.9238795325112867F,  0.9807852804032304F,  0.9951847266721969F,
        0.9987954562051724F,  0.9996988186962042F,  0.9999247018391445F,
        0.9999811752826011F,  0.9999952938095761F,  0.9999988234517018F,
        0.9999997058628822F,  0.9999999264657178F
    };
    static const float Icoef[14] =
    {
        0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
        -0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
        -0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
        -0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
        -0.0007669903187427F, -0.0003834951875714F
    };

    nn = N >> 1;
    ie = N;
    for(n=1; n<=LogN; n++)
    {
        rw = Rcoef[LogN - n];
        iw = Icoef[LogN - n];
        if(Ft_Flag == FT_INVERSE) iw = -iw;
        in = ie >> 1;
        ru = 1.0F;
        iu = 0.0F;
        for(j=0; j<in; j++)
        {
            for(i=j; i<N; i+=ie)
            {
                io       = i + in;
                rtp      = Rdat[i]  + Rdat[io];
                itp      = Idat[i]  + Idat[io];
                rtq      = Rdat[i]  - Rdat[io];
                itq      = Idat[i]  - Idat[io];
                Rdat[io] = rtq * ru - itq * iu;
                Idat[io] = itq * ru + rtq * iu;
                Rdat[i]  = rtp;
                Idat[i]  = itp;
            }

            sr = ru;
            ru = ru * rw - iu * iw;
            iu = iu * rw + sr * iw;
        }

        ie >>= 1;
    }

    for(j=i=1; i<N; i++)
    {
        if(i < j)
        {
            io       = i - 1;
            in       = j - 1;
            rtp      = Rdat[in];
            itp      = Idat[in];
            Rdat[in] = Rdat[io];
            Idat[in] = Idat[io];
            Rdat[io] = rtp;
            Idat[io] = itp;
        }

        k = nn;

        while(k < j)
        {
            j   = j - k;
            k >>= 1;
        }

        j = j + k;
    }

    if(Ft_Flag == FT_DIRECT) return 1;

    rw = 1.0F / N;

    for(i=0; i<N; i++)
    {
        Rdat[i] *= rw;
        Idat[i] *= rw;
    }

    return 1;
}
