// 出来上がったバイナリを実行する時に、関数名をそのまま引数に渡して呼び出すと、テスト用の関数が実行される
// e.g.) ./optimization_test false_sharing_test
//
// all を付けるとすべてのテスト用の関数が実行される
// e.g.) ./optimization_test all

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <mutex>
#include <vector>
#include <memory>
#include <thread>
#include <sstream>
#include <atomic>
#include <chrono>
#include <cstdlib>

struct Timer
{
    using clock_t = std::chrono::high_resolution_clock;
    Timer()
    {
        _begin = clock_t::now();
    }

    double elapsed() const {
        return std::chrono::duration<double>(clock_t::now() - _begin).count();
    }

private:
    clock_t::time_point _begin;
};

void dump(double elapsed_time)
{
    std::printf("Elapsed: %3.8f seconds\n", elapsed_time);
}

template<class T>
struct array2d {
    array2d(int x, int y) : _data(x, std::vector<T>(y)) {}

    std::vector<T> & operator[](int x) { return _data[x]; }
    int size() const { return _data.size(); }

    bool operator==(array2d const & rhs) {
        if(size() != rhs.size()) {
            return false;
        }

        for(int i = 0, end = size(); i < end; ++i) {
            if(_data[i] != rhs._data[i]) {
                return false;
            }
        }

        return true;
    }

    bool operator!=(array2d const & rhs) {
        return !(*this == rhs);
    }
private:
    std::vector<std::vector<T>> _data;
};

template<class T>
struct array3d {
    array3d(int x, int y, int z) : _data(x, array2d<T>(y, z)) {}

    array2d<T> & operator[](int x) { return _data[x]; }

    int size() const { return _data.size(); }

    bool operator==(array3d const & rhs) {
        if(size() != rhs.size()) {
            return false;
        }

        for(int i = 0, end = size(); i < end; ++i) {
            if(_data[i] != rhs._data[i]) {
                return false;
            }
        }

        return true;
    }

    bool operator!=(array3d const & rhs) {
        return !(*this == rhs);
    }

private:
    std::vector<array2d<T>> _data;
};

void access_order_test();
void use_loop_blocking_test();
void use_loop_blocking_test_for_sample_interleaving();
void cache_thrashing_test();
void false_sharing_test();
void abs_impl_test();

int main(int argc, char **argv)
{
    if(argc != 2) { return 0; }

    auto test_case_name = std::string(argv[1]);

    bool const all = (test_case_name == "all");

    if(all || test_case_name == "access_order_test") {
        access_order_test();
    }
    if(all || test_case_name == "use_loop_blocking_test") {
        use_loop_blocking_test();
    }
    if(all || test_case_name == "use_loop_blocking_test_for_sample_interleaving") {
        use_loop_blocking_test_for_sample_interleaving();
    }
    if(all || test_case_name == "cache_thrashing_test") {
        cache_thrashing_test();
    }
    if(all || test_case_name == "false_sharing_test") {
        false_sharing_test();
    }
    if(all || test_case_name == "abs_impl_test") {
        abs_impl_test();
    }
}

// メモリアクセスは、飛び飛びにアクセスするよりも、
// 連続した領域に順番にアクセスしていくほうが速いというのを検証するテスト
// このテストでは行列を表す多次元配列にアクセスする際に、どの次元で要素をイテレートするかを変えて
// 負荷を計測する。
void access_order_test()
{
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "access_order_test" << std::endl;
    std::cout << std::endl;

    // http://cocodrips.hateblo.jp/entry/2014/01/26/134501 を改変

    int const kSize = 500;
    array3d<int> arr(kSize, kSize, kSize);

    int sum = 0;
    std::cout << "Preferable Order" << std::endl;
    {
        sum = 0;
        Timer t;
        for (int i = 0; i < kSize; i++) {
            for (int j = 0; j < kSize; j++) {
                for (int k = 0; k < kSize; k++) {
                    // arr の末尾の配列はメモリ上に連続的に並んでいるので、
                    // これに対して連続的にアクセスするのは効率的
                    arr[i][j][k] = i * j * k;
                    sum += arr[i][j][k];
                }
            }
        }
        dump(t.elapsed());
    }
    std::cout << "sum : " << sum << std::endl;

    std::cout << "Wrong Order" << std::endl;
    {
        sum = 0;
        Timer t;
        for (int k = 0; k < kSize; k++) {
            for (int j = 0; j < kSize; j++) {
                for (int i = 0; i < kSize; i++) {
                    // こちらのコードでは、ある行列のi番目の要素にある配列を取得し、
                    // その中のj番目の要素を取得し、その中のk番目の要素を取得するというのを
                    // 毎回i番目の要素の取得から行う。
                    // そのため効率が悪い。
                    arr[i][j][k] = i * j * k;
                    sum += arr[i][j][k];
                }
            }
        }
        dump(t.elapsed());
    }
    std::cout << "sum : " << sum << std::endl;
}

// ループ処理を行う際に、処理単位を分割することでキャッシュ効率を向上させる手法のテスト
// 行列の積の計算は、行アクセスと列アクセスが混在するので、片方を効率よく処理させようとすると
// もう片方のアクセス効率が落ちる問題がある。
// これに対して、積の計算をするときの処理単位を小さく分割することでそれぞれの効率を下げずに計算できる。
void use_loop_blocking_test()
{
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "use_loop_blocking_test" << std::endl;
    std::cout << std::endl;

    int const kSize = 150;

    array2d<std::int64_t> a(kSize, kSize);
    array2d<std::int64_t> b(kSize, kSize);
    array2d<std::int64_t> c(kSize, kSize);
    array2d<std::int64_t> d(kSize, kSize);

    for(int i = 0; i < kSize; ++i) {
        for(int j = 0; j < kSize; ++j) {
            a[i][j] = i + j;
            b[i][j] = i + j + 2;
        }
    }

    std::cout << "Don't use loop blocking" << std::endl;
    {
        Timer t;
        for(int i = 0; i < kSize; ++i) {
            for(int j = 0; j < kSize; ++j) {
                c[i][j] = 0;
            }
        }

        for (int i = 0; i < kSize; i++) {
            for (int j = 0; j < kSize; j++) {
                for (int k = 0; k < kSize; ++k) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        dump(t.elapsed());
    }

    std::cout << "Use loop blocking" << std::endl;
    {
        int fastest_block_size = 1000000;
        double fastest_block_time = 1000000;

        for(int block_size = 10; block_size < 256; ++block_size) {
            std::cout << "block_size: " << block_size << std::endl;

            Timer t;
            for(int i = 0; i < kSize; ++i) {
                for(int j = 0; j < kSize; ++j) {
                    d[i][j] = 0;
                }
            }

            {
                for(int ib = 0; ib < kSize; ib += block_size) {
                    int const iend = std::min<int>(ib + block_size, kSize);
                    for(int jb = 0; jb < kSize; jb += block_size) {
                        int const jend = std::min<int>(jb + block_size, kSize);
                        for(int kb = 0; kb < kSize; kb += block_size) {
                            int const kend = std::min<int>(kb + block_size, kSize);
                            for (int i = 0; i < iend; i++) {
                                for (int j = 0; j < jend; j++) {
                                    for (int k = 0; k < kend; ++k) {
                                        d[i][j] += a[i][k] * b[k][j];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            auto const e = t.elapsed();
            dump(e);

            if(fastest_block_time > e) {
                fastest_block_time = e;
                fastest_block_size = block_size;
            }

            if(c == d) {
                std::cout << "incorrect result." << std::endl;
            }
        }

        std::cout << std::endl;
        std::cout << "fastest block size: " << fastest_block_size << std::endl;
        dump(fastest_block_time);
    }
}

// wavファイルのデータ形式のように、1サンプルデータをチャンネル数で並べ、それをサンプル時間分並べたデータ構造と、
// プログラム上で扱いやすいような、チャンネル数の配列の中に、それぞれサンプル時間分のサンプルが並ぶデータ構造では、
// メモリアクセスする順番がことなるため、どちらかだけで連続的に処理を行うと、もう一方は常にキャッシュ効率が悪くなってしまう
// そこで処理する単位を分割し、どちらもキャッシュ効率を高めるようにすると、パフォーマンスが向上すると考えられる。
// それを検証するテスト
void use_loop_blocking_test_for_sample_interleaving()
{
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "use_loop_blocking_test_for_sample_interleaving" << std::endl;
    std::cout << std::endl;

    int const kSize = 96000 * 60 * 5;
    array2d<std::int16_t> d(kSize, 2); // wavfile
    array2d<float> e(2, kSize); // プログラム上で持ちたい形式
    int sum;

    std::cout << "Don't use loop blocking" << std::endl;
    {
        sum = 0;
        Timer t;
        for (int smp = 0; smp < kSize; smp++) {
            for (int ch = 0; ch < 2; ch++) {
                e[ch][smp] = d[smp][ch] / 32768.0;
                sum += e[ch][smp];
            }
        }
        dump(t.elapsed());
    }

    std::cout << "Estimate preferable blocking size for sample interleaving." << std::endl;
    {
        int fastest_block_size = 1000000;
        double fastest_block_time = 1000000;

        for(int block_size = 1; block_size < 128; ++block_size) {
            std::cout << "block_size: " << block_size << std::endl;
            sum = 0;
            Timer t;
            {
                for (int smpb = 0; smpb < kSize; smpb += block_size) {
                    int smpend = std::min<int>(smpb + block_size, kSize);
                    for(int smp = smpb; smp < smpend; ++smp) {
                        for(int ch = 0; ch < 2; ++ch) {
                            e[ch][smp] = d[smp][ch] / 32768.0;
                            sum += e[ch][smp];
                        }
                    }
                }
            }

            auto const e = t.elapsed();
            dump(e);

            if(fastest_block_time > e) {
                fastest_block_time = e;
                fastest_block_size = block_size;
            }
        }

        std::cout << std::endl;
        std::cout << "fastest block size: " << fastest_block_size << std::endl;
        dump(fastest_block_time);
    }

}

namespace cache_thrashing_test_data {
    int const kSize = 2048;
    int const kWay = 16;

    struct MayOccur {
        double a[kSize];
        double b[kSize];
    } mo[kWay];

    struct MayNotOccur
    {
        double c[kSize];
        double dummy[64];
        double d[kSize];
    } mno[kWay];
}

// CPUがデータをキャッシュする単位は、バイト単位ではなく、キャッシュラインと呼ばれる一定の大きさのデータのかたまりである
// CPUが、メモリ上のある変数Xの値を読み込み(1)、そこから2のべき乗だけ離れた位置にある変数Yに値を書き込む(2)と以下のような問題が起きる
// (1)のタイミングでXを含むメモリ領域がキャッシュに載せられたあと、
// (2)のタイミングでYを含むメモリ領域がキャッシュに載せられようとするが、
// 2のべき乗だけ離れた位置にあるメモリは、キャッシュラインの位置がかぶってしまう（セットアソシエイティブ方式）
// そのためCPUはXを含むメモリ領域をキャッシュから退避し、Yを新たにキャッシュに載せる
// その後また(1)が発生すると、CPUはYを含むキャッシュを退避し、Xを新たにキャッシュに載せる
// このような状態が続くと、毎回キャッシュが無効化され、キャッシュ効率が低下してしまう。
// これをキャッシュスラッシングと呼ぶ
// これにはCPU側で対策がされていて、キャッシュラインのテーブルを複数持つようにし、
// XとYのアドレスからキャッシュラインの位置を計算したときに同じ位置になってしまっても、
// それぞれを別のテーブルに割り当てるようにし、お互いのキャッシュラインが無効化されることがないようになっている
//（Nウェイセットアソシエイティブ方式）
// ただし、このテーブル数を超えたサイズのキャッシュがキャッシュラインに載せられようとしたときには、いずれかのキャッシュラインが無効化されるので、キャッシュ効率が低下する可能性がある。
// これにプログラム側で対処するには、大量にアクセスするいくつかのメモリ領域が、それぞれ2のべき乗倍だけ離れた位置にならないようにする。
// このプログラムでは、変数c, dの間にダミーの変数を設けて、アクセスする位置の差が2のべき乗倍にならないようにしている
void cache_thrashing_test()
{
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "cache_thrashing_test" << std::endl;
    std::cout << std::endl;

    using namespace cache_thrashing_test_data;

    int const kRepeat = 50000;
    int sum = 0;

    std::cout << "Cache thrashing may not occur" << std::endl;
    {
        sum = 0;
        Timer t;
        for(int rep = 0; rep < kRepeat; ++rep) {
            for(int i = 0; i < kSize; ++i) {
                for(int w = 0; w < kWay; ++w) {
                    mno[w].c[i] += mno[w].d[i];
                    sum += mno[w].c[i];
                }
            }
        }
        dump(t.elapsed());
    }

    std::cout << "Cache thrashing may occur" << std::endl;
    {
        sum = 0;
        Timer t;
        for(int rep = 0; rep < kRepeat; ++rep) {
            for(int i = 0; i < kSize; ++i) {
                for(int w = 0; w < kWay; ++w) {
                    mo[w].a[i] += mo[w].b[i];
                    sum += mo[w].a[i];
                }
            }
        }
        dump(t.elapsed());
    }
}

namespace false_sharing_test_data
{
    int const kNumThreads = 16;
    int const kSkip = 16;

    struct FalseSharingMayNotOccur
    {
        int a_[kNumThreads * kSkip];
    } mno;

    struct FalseSharingMayOccur
    {
        int a_[kNumThreads];
    } mo;
}

// あるCPUコア1が変数Xにアクセスし(1)、その値を変更すると(2)、
// (1)のタイミングでXを含む周辺がコア1のキャッシュラインに載せられ、
// (2)のタイミングで、そのキャッシュラインのデータが変更される。
// このとき、他のコア2もXの近傍にある別の変数Yにアクセスしていたとしよう。
// そのとき、Yの周辺（そしてXも含まれる、(1)と同じメモリ領域）がコア2のキャッシュラインに載せられていることになる。
// (2)のときに、コア1でキャッシュライン内のデータが変更されると、
// CPUは他のコアに対してそのキャッシュラインが変更されたことを通知し、コア2のキャッシュラインは無効化される。
// そのため、コア2でYに再びアクセスしようとすると、Yの周辺を再びコア2のキャッシュラインに載せなければならなくなる。
// 逆にコア2でYの値を変更すると、コア1のキャッシュラインが無効化される
// これを繰り返すと、スレッド間でお互いにキャッシュラインを無効化し合うことになり、キャッシュ効率が低下してしまう。
// この問題を、フォールスシェアリング（偽の共有）と呼ぶ。
// この問題の面白い所は、ソース上ではコア1とコア2は別の位置にあるオブジェクトへアクセスしているはずなのに、それらが互いに影響を与えてしまうところである。
// この問題を防ぐには、コア1とコア2でアクセスするそれぞれのデータが、キャッシュラインのサイズより遠くなるようにする
// このコードではkSkipという変数でそれぞれのコアがアクセスするオブジェクトの位置を離すことで、フォールスシェアリングを回避している
void false_sharing_test()
{
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "false_sharing_test" << std::endl;
    std::cout << std::endl;

    using namespace false_sharing_test_data;

    int const kRepeat = 4000;
    int sum1 = 0;
    int sum2 = 0;

    std::cout << "False sharing may not occur" << std::endl;
    {
        std::vector<std::thread> ths(kNumThreads);

        std::mutex mtx;
        double time_sum = 0;

        for(int n = 0; n < kNumThreads; ++n) {
            ths[n] = std::thread([&, n] {
                Timer t;
                for(int i = 0; i < kRepeat; ++i) {
                    for(int j = 0; j < kRepeat; ++j) {
                        mno.a_[n * kSkip] += (i * j);
                        sum1 += mno.a_[n * kSkip];
                    }
                }

                auto const e = t.elapsed();

                std::lock_guard<std::mutex> lock(mtx);
                time_sum += e;
            });
        }

        for(int n = 0; n < kNumThreads; ++n) {
            ths[n].join();
        }

        dump(time_sum / kNumThreads);
    }

    std::cout << "False sharing may occur" << std::endl;
    {
        std::vector<std::thread> ths(kNumThreads);

        std::mutex mtx;
        double time_sum = 0;

        for(int n = 0; n < kNumThreads; ++n) {
            ths[n] = std::thread([&, n] {
                Timer t;
                for(int i = 0; i < kRepeat; ++i) {
                    for(int j = 0; j < kRepeat; ++j) {
                        mo.a_[n] += (i * j);
                        sum1 += mo.a_[n];
                    }
                }

                auto const e = t.elapsed();

                std::lock_guard<std::mutex> lock(mtx);
                time_sum += e;
            });
        }

        for(int n = 0; n < kNumThreads; ++n) {
            ths[n].join();
        }

        dump(time_sum / kNumThreads);
    }
}

int64_t abs_using_if(int64_t x)
{
    return x >= 0 ? x : -x;
};

int64_t abs_using_bitop(int64_t x)
{
    std::uint64_t M = (x >> 63);
    return (x ^ M) - M;
};

// abs()を実装するとき、条件分岐を使って実装するよりもビット演算だけで処理を完結させるほうが速いように思われる。
// しかし実際には、どちらもさほど変わらない。
// Intelの最近の命令セットには、フラグの状態に合わせてmovを行う命令があるので、条件分岐 + movを1命令で行える。
// したがって、最適化によってどちらも同じような性能になる
void abs_impl_test()
{
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "fast_abs_test" << std::endl;
    std::cout << std::endl;

    assert(abs_using_if(-(0xF1F2F3F4F5LL)) == abs(0xF1F2F3F4F5LL));
    assert(abs_using_bitop(-(0xF1F2F3F4F5LL)) == abs(0xF1F2F3F4F5LL));

    int64_t count = 0;
    {
        // 条件分岐を伴う abs の実装

        std::cout << "get absolute value using if" << std::endl;
        count = 0;
        Timer t;
        {
            for(int64_t i = 0; i < 1000 * 1000 * 100; ++i) {
                int64_t tmp = -(i & 1) * i + ((i+1) & 1) * i;
                count += (count & i) + abs_using_if(tmp);
            }
        }
        dump(t.elapsed());
        std::cout << count << std::endl;
    }
    {
        // bit 演算による abs の実装

        std::cout << "get absolute value using bit operation" << std::endl;
        count = 0;
        Timer t;
        {
            for(int64_t i = 0; i < 1000 * 1000 * 100; ++i) {
                int64_t tmp = -(i & 1) * i + ((i+1) & 1) * i;
                count += (count & i) + abs_using_bitop(tmp);
            }
        }
        dump(t.elapsed());
        std::cout << count << std::endl;
    }
}
