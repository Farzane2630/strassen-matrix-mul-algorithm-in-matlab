# Strassen's Matrix Multiplication Algorithm

### Algorithm Design & Analysis — University Coursework

**Implementation:** MATLAB / GNU Octave  
**Algorithm:** Strassen's Recursive Matrix Multiplication  
**Purpose:** Understanding, implementing, and comparing Strassen's algorithm with conventional matrix multiplication.

---

# 🇬🇧

## 1. Project Overview

This project implements **Strassen's Matrix Multiplication Algorithm**, a divide-and-conquer algorithm that multiplies two square matrices using fewer recursive multiplications than the conventional matrix multiplication method.

The main objectives are:

- Understand the mathematical foundation of Strassen's algorithm.
- Learn how divide-and-conquer works in matrix multiplication.
- Implement the algorithm using MATLAB / GNU Octave.
- Compare its results with MATLAB's built-in matrix multiplication.
- Measure and analyze execution time.
- Investigate how matrix size and the recursion cutoff affect performance.

## 2. What Is Matrix Multiplication?

Matrix multiplication is a mathematical operation that combines two matrices to produce a third matrix.

Suppose we have two square matrices:

\[
A =
\begin{bmatrix}
a*{11} & a*{12}\\
a*{21} & a*{22}
\end{bmatrix}
\qquad
B =
\begin{bmatrix}
b*{11} & b*{12}\\
b*{21} & b*{22}
\end{bmatrix}
\]

Their product is:

\[
C = A \times B
\]

Each element of the resulting matrix is calculated by multiplying the corresponding row of A by the corresponding column of B and adding the results.

For example:

\[
C*{11}=a*{11}b*{11}+a*{12}b\_{21}
\]

The complete result is:

\[
C =
\begin{bmatrix}
a*{11}b*{11}+a*{12}b*{21} &
a*{11}b*{12}+a*{12}b*{22}\\
a*{21}b*{11}+a*{22}b*{21} &
a*{21}b*{12}+a*{22}b*{22}
\end{bmatrix}
\]

For two \(n \times n\) matrices, conventional matrix multiplication requires approximately:

\[
T(n)=O(n^3)
\]

This is because the standard algorithm performs \(n^3\) scalar multiplications and additions at the same asymptotic order.

## 3. What Is Strassen's Algorithm?

Strassen's algorithm was introduced by German mathematician Volker Strassen in 1969.

It was an important discovery because it demonstrated that matrix multiplication could be performed with fewer than eight recursive multiplications when two matrices are divided into four submatrices.

The key idea is:

**Instead of performing eight smaller matrix multiplications, Strassen's algorithm performs only seven.**

This reduction improves the theoretical time complexity.

### 3.1 Divide-and-Conquer Approach

Strassen's algorithm uses the divide-and-conquer strategy.

The general process consists of three steps:

1. **Divide:** Split each input matrix into four smaller submatrices.
2. **Conquer:** Recursively multiply the submatrices using Strassen's formulas.
3. **Combine:** Use the seven intermediate products to construct the final result.

Suppose A and B are \(n \times n\) matrices:

\[
A =
\begin{bmatrix}
A*{11} & A*{12}\\
A*{21} & A*{22}
\end{bmatrix}
\qquad
B =
\begin{bmatrix}
B*{11} & B*{12}\\
B*{21} & B*{22}
\end{bmatrix}
\]

Each submatrix has dimensions:

\[
\frac{n}{2}\times\frac{n}{2}
\]

The final result is also divided into four blocks:

\[
C =
\begin{bmatrix}
C*{11} & C*{12}\\
C*{21} & C*{22}
\end{bmatrix}
\]

### 3.2 The Seven Multiplications

This is the most important mathematical part of Strassen's algorithm.

Instead of calculating eight products, we calculate seven intermediate matrices:

\[
\begin{aligned}
M*1 &= (A*{11}+A*{22})(B*{11}+B*{22})\\
M_2 &= (A*{21}+A*{22})B*{11}\\
M*3 &= A*{11}(B*{12}-B*{22})\\
M*4 &= A*{22}(B*{21}-B*{11})\\
M*5 &= (A*{11}+A*{12})B*{22}\\
M*6 &= (A*{21}-A*{11})(B*{11}+B*{12})\\
M_7 &= (A*{12}-A*{22})(B*{21}+B\_{22})
\end{aligned}
\]

These seven products are then combined to calculate the four blocks of the resulting matrix:

\[
C\_{11}=M_1+M_4-M_5+M_7
\]

\[
C\_{12}=M_3+M_5
\]

\[
C\_{21}=M_2+M_4
\]

\[
C\_{22}=M_1-M_2+M_3+M_6
\]

Finally:

\[
C =
\begin{bmatrix}
C*{11} & C*{12}\\
C*{21} & C*{22}
\end{bmatrix}
\]

**Why is this faster?**

Conventional divide-and-conquer multiplication requires eight recursive multiplications. Strassen requires only seven.

Although Strassen performs additional matrix additions and subtractions, reducing the number of recursive multiplications improves the asymptotic complexity.

## 4. Time Complexity Analysis

The conventional matrix multiplication algorithm has the following complexity:

\[
T(n)=O(n^3)
\]

For Strassen's algorithm, each recursive call performs:

- 7 multiplications of matrices with size \(n/2\).
- A number of matrix additions and subtractions with total cost \(O(n^2)\).

Its recurrence relation is:

\[
T(n)=7T(n/2)+O(n^2)
\]

Using the Master Theorem:

\[
T(n)=O(n^{\log_2 7})
\]

Since:

\[
\log_2 7 \approx 2.807
\]

The final time complexity is:

\[
\boxed{T(n)=O(n^{2.807})}
\]

Comparison:

| Algorithm                   | Time Complexity  |
| --------------------------- | ---------------- |
| Conventional multiplication | \(O(n^3)\)       |
| Strassen's algorithm        | \(O(n^{2.807})\) |

Strassen's algorithm has a better theoretical asymptotic complexity. However, this does not mean it is always faster in practice. Recursion, memory allocation, and matrix additions introduce overhead.

## 5. Project Structure

The project contains two MATLAB files:

```text
Algo-design-strassen-algo/
│
├── strassen_demo.m
│
└── strassen_multiply.m
```

### File 1: `strassen_multiply.m`

This file contains the main Strassen multiplication implementation.

Its responsibilities include:

- Receiving two input matrices and a cutoff value.
- Checking matrix dimensions.
- Padding matrices with zeros when necessary.
- Calling the recursive Strassen implementation.
- Returning the final multiplication result.

**Padding:**

Strassen's recursive implementation is easiest to apply when the matrix dimensions are powers of two.

For example, if the input matrices are \(10 \times 10\), the implementation can extend them to \(16 \times 16\) by adding zeros.

The padded matrices look conceptually like this:

\[
A' =
\begin{bmatrix}
A & 0\\
0 & 0
\end{bmatrix}
\]

The additional rows and columns contain zeros, so the original multiplication result can be recovered by extracting the original \(10 \times 10\) portion.

### File 2: `strassen_demo.m`

This file is the demonstration and testing script.

It performs the following operations:

1. Clears the Command Window and workspace.
2. Defines the matrix size and recursion cutoff.
3. Generates two random matrices.
4. Displays the input matrices.
5. Executes Strassen's algorithm.
6. Measures its execution time.
7. Performs conventional multiplication using MATLAB's built-in `*` operator.
8. Measures the built-in multiplication time.
9. Displays both results.
10. Calculates the maximum absolute error between the two results.

## 6. Understanding the Cutoff Parameter

The implementation uses a parameter called `cutoff`.

For example:

```matlab
n = 10;
cutoff = 2;
```

Strassen's algorithm uses recursion to divide matrices into smaller submatrices.

However, continuing recursion for very small matrices is inefficient because recursive function calls and matrix operations introduce additional overhead.

The cutoff determines when the algorithm stops dividing matrices and switches to conventional multiplication.

For example:

```matlab
cutoff = 2;
```

When the recursive matrix size reaches the cutoff threshold, ordinary multiplication is used as the base case.

This technique is called a **hybrid approach** because it combines Strassen's algorithm with conventional matrix multiplication.

Experimenting with different cutoff values can help identify which configuration performs best for a particular computer and matrix size.

## 7. Running the Project

### Requirements

- MATLAB or GNU Octave.
- Both `.m` files in the same directory.

No external libraries are required.

### Using GNU Octave

Open the terminal and navigate to the project directory:

```bash
cd ~/Documents/UNI-semester-3/courseworks/Algo-design-strassen-algo
```

Start Octave:

```bash
octave
```

Execute the demonstration:

```matlab
strassen_demo
```

Alternatively, execute the script directly from the Linux terminal:

```bash
octave --no-gui --quiet strassen_demo.m
```

### Using MATLAB

Open MATLAB, navigate to the project folder, and execute:

```matlab
strassen_demo
```

**Important:** Both files must be accessible from MATLAB / Octave. The function file must be named `strassen_multiply.m` because the main script calls that function by name.

## 8. Performance Comparison

The demonstration compares two approaches:

**Approach 1: Strassen's algorithm**

```matlab
tic;
C_strassen = strassen_multiply(A, B, cutoff);
t_strassen = toc;
```

**Approach 2: Built-in matrix multiplication**

```matlab
tic;
C_builtin = A * B;
t_builtin = toc;
```

The `tic` and `toc` commands measure execution time.

The program reports:

- Strassen execution time.
- Built-in multiplication execution time.
- Maximum absolute error.

### Why might built-in multiplication be faster?

For small matrices, MATLAB's built-in multiplication may be faster because it uses highly optimized numerical libraries and avoids the overhead of recursive function calls.

Strassen's theoretical advantage becomes more relevant for sufficiently large matrices.

Therefore, a small matrix experiment does not necessarily demonstrate Strassen's performance advantage.

For more meaningful experiments, try different matrix sizes, such as:

- 10 × 10
- 100 × 100
- 500 × 500
- 1000 × 1000

Also experiment with different cutoff values.

For reliable timing comparisons, repeat each experiment several times and compare average execution times.

## 9. Numerical Accuracy

The program calculates the maximum absolute error using:

```matlab
max(abs(C_strassen(:) - C_builtin(:)))
```

The colon operator `(:)` converts a matrix into a column vector.

The `abs()` function calculates the absolute difference between corresponding elements.

The `max()` function finds the largest difference.

Mathematically:

\[
E*{\max} =
\max*{i,j}
|C*{\text{Strassen}}(i,j)-C*{\text{Builtin}}(i,j)|
\]

If the maximum error is very close to zero, the two implementations produce numerically equivalent results.

For example:

```text
Maximum Error: 2.22e-15
```

This is a very small floating-point difference and is generally expected in numerical computations.

It does not necessarily indicate a mistake in the algorithm.

## 10. Expected Output

The demonstration displays:

```text
========== FINAL RESULT ==========

Strassen Result:
[Matrix C]

Built-in Result:
[Matrix C]

Strassen Time: ...
Built-in Time: ...
Maximum Error: ...
```

The actual execution times and numerical error depend on the computer, software environment, matrix size, and cutoff value.

## 11. Key Takeaways

After studying and running this project, students should understand:

1. How conventional matrix multiplication works.
2. How divide-and-conquer can be applied to matrix multiplication.
3. Why Strassen's algorithm uses seven recursive multiplications instead of eight.
4. How seven intermediate products produce the final four matrix blocks.
5. How recursion and base cases work.
6. Why padding is useful for matrix dimensions that are not powers of two.
7. Why the theoretical complexity is \(O(n^{2.807})\).
8. Why a better theoretical complexity does not guarantee better performance for every input size.
9. How to compare two algorithms experimentally.
10. How floating-point precision affects numerical results.

---

# 🇮🇷

## ۱. معرفی پروژه

این پروژه پیاده‌سازی **الگوریتم ضرب ماتریس استراسن (Strassen's Matrix Multiplication Algorithm)** است.

الگوریتم استراسن یک الگوریتم مبتنی بر روش تقسیم و حل (Divide and Conquer) است که با استفاده از تعداد کمتری ضرب بازگشتی نسبت به روش معمول ضرب ماتریس‌ها، دو ماتریس مربعی را در یکدیگر ضرب می‌کند.

اهداف اصلی این پروژه عبارت‌اند از:

- درک مبانی ریاضی الگوریتم استراسن.
- آشنایی با روش تقسیم و حل در ضرب ماتریس‌ها.
- پیاده‌سازی الگوریتم با استفاده از MATLAB یا GNU Octave.
- مقایسه نتیجه با عملگر ضرب داخلی MATLAB.
- اندازه‌گیری و تحلیل زمان اجرا.
- بررسی تأثیر اندازه ماتریس و پارامتر cutoff بر عملکرد الگوریتم.

## ۲. ضرب ماتریس چیست؟

ضرب ماتریس یک عملیات ریاضی است که دو ماتریس را با یکدیگر ترکیب می‌کند و یک ماتریس جدید به وجود می‌آورد.

فرض کنید دو ماتریس مربعی زیر را داریم:

\[
A =
\begin{bmatrix}
a*{11} & a*{12}\\
a*{21} & a*{22}
\end{bmatrix}
\qquad
B =
\begin{bmatrix}
b*{11} & b*{12}\\
b*{21} & b*{22}
\end{bmatrix}
\]

حاصل‌ضرب آن‌ها به شکل زیر محاسبه می‌شود:

\[
C = A \times B
\]

برای محاسبه هر عنصر از ماتریس حاصل، باید سطر متناظر از ماتریس اول را در ستون متناظر از ماتریس دوم ضرب داخلی کنیم.

برای مثال:

\[
C*{11}=a*{11}b*{11}+a*{12}b\_{21}
\]

در نتیجه:

\[
C =
\begin{bmatrix}
a*{11}b*{11}+a*{12}b*{21} &
a*{11}b*{12}+a*{12}b*{22}\\
a*{21}b*{11}+a*{22}b*{21} &
a*{21}b*{12}+a*{22}b*{22}
\end{bmatrix}
\]

در روش معمول ضرب دو ماتریس با ابعاد \(n \times n\)، پیچیدگی زمانی برابر است با:

\[
T(n)=O(n^3)
\]

زیرا این روش به ترتیب \(n^3\) ضرب عددی و تعداد مشابهی عملیات جمع نیاز دارد.

## ۳. الگوریتم استراسن چیست؟

الگوریتم استراسن در سال ۱۹۶۹ توسط ریاضی‌دان آلمانی، فولکر استراسن (Volker Strassen)، معرفی شد.

این الگوریتم یک کشف مهم در علم طراحی الگوریتم بود؛ زیرا نشان داد که برای ضرب ماتریس‌ها می‌توان از تعداد ضرب بازگشتی کمتری نسبت به روش تقسیم و حل معمول استفاده کرد.

ایده اصلی الگوریتم این است:

**به‌جای انجام هشت ضرب بازگشتی روی ماتریس‌های کوچک‌تر، تنها هفت ضرب انجام می‌دهیم.**

این کاهش تعداد ضرب‌ها باعث بهبود پیچیدگی زمانی الگوریتم می‌شود.

### ۳.۱. روش تقسیم و حل (Divide and Conquer)

الگوریتم استراسن از روش تقسیم و حل استفاده می‌کند.

این روش شامل سه مرحله اصلی است:

۱. **تقسیم (Divide):** هر ماتریس ورودی به چهار زیرماتریس کوچک‌تر تقسیم می‌شود.

۲. **حل (Conquer):** ضرب زیرماتریس‌ها با استفاده از فرمول‌های استراسن و به‌صورت بازگشتی انجام می‌شود.

۳. **ترکیب (Combine):** هفت حاصل‌ضرب میانی با یکدیگر ترکیب می‌شوند تا ماتریس نهایی ساخته شود.

فرض کنید دو ماتریس A و B با ابعاد \(n \times n\) داریم:

\[
A =
\begin{bmatrix}
A*{11} & A*{12}\\
A*{21} & A*{22}
\end{bmatrix}
\qquad
B =
\begin{bmatrix}
B*{11} & B*{12}\\
B*{21} & B*{22}
\end{bmatrix}
\]

هر زیرماتریس دارای ابعاد زیر است:

\[
\frac{n}{2}\times\frac{n}{2}
\]

ماتریس حاصل نیز به چهار بلوک تقسیم می‌شود:

\[
C =
\begin{bmatrix}
C*{11} & C*{12}\\
C*{21} & C*{22}
\end{bmatrix}
\]

### ۳.۲. هفت ضرب اصلی الگوریتم استراسن

این بخش مهم‌ترین قسمت ریاضی الگوریتم استراسن محسوب می‌شود.

به‌جای محاسبه هشت ضرب معمول، هفت ماتریس میانی را محاسبه می‌کنیم:

\[
\begin{aligned}
M*1 &= (A*{11}+A*{22})(B*{11}+B*{22})\\
M_2 &= (A*{21}+A*{22})B*{11}\\
M*3 &= A*{11}(B*{12}-B*{22})\\
M*4 &= A*{22}(B*{21}-B*{11})\\
M*5 &= (A*{11}+A*{12})B*{22}\\
M*6 &= (A*{21}-A*{11})(B*{11}+B*{12})\\
M_7 &= (A*{12}-A*{22})(B*{21}+B\_{22})
\end{aligned}
\]

سپس از این هفت حاصل‌ضرب برای محاسبه چهار بلوک ماتریس نهایی استفاده می‌کنیم:

\[
C\_{11}=M_1+M_4-M_5+M_7
\]

\[
C\_{12}=M_3+M_5
\]

\[
C\_{21}=M_2+M_4
\]

\[
C\_{22}=M_1-M_2+M_3+M_6
\]

در نهایت ماتریس حاصل به شکل زیر ساخته می‌شود:

\[
C =
\begin{bmatrix}
C*{11} & C*{12}\\
C*{21} & C*{22}
\end{bmatrix}
\]

**چرا این روش سریع‌تر است؟**

در روش تقسیم و حل معمول برای ضرب ماتریس‌ها، هشت ضرب بازگشتی انجام می‌شود؛ اما استراسن تنها به هفت ضرب نیاز دارد.

اگرچه در الگوریتم استراسن تعداد عملیات جمع و تفریق ماتریسی افزایش پیدا می‌کند، کاهش تعداد ضرب‌های بازگشتی باعث بهبود پیچیدگی زمانی می‌شود.

## ۴. تحلیل پیچیدگی زمانی (Time Complexity)

پیچیدگی زمانی الگوریتم معمول ضرب ماتریس‌ها برابر است با:

\[
T(n)=O(n^3)
\]

اما در الگوریتم استراسن، در هر مرحله:

- هفت ضرب بازگشتی روی ماتریس‌هایی با اندازه \(n/2\) انجام می‌شود.
- تعدادی عملیات جمع و تفریق ماتریسی با مجموع هزینه \(O(n^2)\) انجام می‌شود.

رابطه بازگشتی الگوریتم به شکل زیر است:

\[
T(n)=7T(n/2)+O(n^2)
\]

با استفاده از قضیه مستر (Master Theorem):

\[
T(n)=O(n^{\log_2 7})
\]

از آنجا که:

\[
\log_2 7 \approx 2.807
\]

پیچیدگی زمانی نهایی برابر است با:

\[
\boxed{T(n)=O(n^{2.807})}
\]

مقایسه دو الگوریتم:

| الگوریتم            | پیچیدگی زمانی    |
| ------------------- | ---------------- |
| ضرب معمول ماتریس‌ها | \(O(n^3)\)       |
| الگوریتم استراسن    | \(O(n^{2.807})\) |

الگوریتم استراسن از نظر پیچیدگی زمانی تئوری عملکرد بهتری دارد؛ اما این به معنی سریع‌تر بودن آن در تمام شرایط نیست. فراخوانی‌های بازگشتی، تخصیص حافظه و عملیات جمع و تفریق می‌توانند سربار محاسباتی ایجاد کنند.

## ۵. ساختار پروژه

این پروژه شامل دو فایل MATLAB است:

```text
Algo-design-strassen-algo/
│
├── strassen_demo.m
│
└── strassen_multiply.m
```

### فایل اول: `strassen_multiply.m`

این فایل شامل پیاده‌سازی اصلی الگوریتم ضرب ماتریس استراسن است.

وظایف این فایل عبارت‌اند از:

- دریافت دو ماتریس ورودی و مقدار cutoff.
- بررسی ابعاد ماتریس‌ها.
- اضافه کردن صفر به ماتریس‌ها در صورت نیاز.
- اجرای پیاده‌سازی بازگشتی الگوریتم استراسن.
- بازگرداندن نتیجه نهایی ضرب ماتریس‌ها.

**Padding یا اضافه کردن صفر:**

پیاده‌سازی بازگشتی استراسن زمانی ساده‌تر است که ابعاد ماتریس‌ها توان‌های عدد ۲ باشند.

برای مثال، اگر ماتریس ورودی ابعاد \(10 \times 10\) داشته باشد، برنامه می‌تواند با اضافه کردن سطرها و ستون‌های صفر، ابعاد آن را به \(16 \times 16\) تبدیل کند.

ماتریس جدید به‌صورت مفهومی چنین خواهد بود:

\[
A' =
\begin{bmatrix}
A & 0\\
0 & 0
\end{bmatrix}
\]

از آنجا که سطرها و ستون‌های اضافه‌شده شامل صفر هستند، نتیجه ضرب ماتریس اصلی با استخراج قسمت اولیه ماتریس حاصل قابل بازیابی است.

### فایل دوم: `strassen_demo.m`

این فایل اسکریپت اصلی برای اجرا و آزمایش الگوریتم است.

وظایف آن عبارت‌اند از:

۱. پاک کردن محیط Command Window و متغیرهای قبلی.

۲. تعیین اندازه ماتریس‌ها و مقدار cutoff.

۳. تولید دو ماتریس تصادفی.

۴. نمایش ماتریس‌های ورودی.

۵. اجرای الگوریتم استراسن.

۶. اندازه‌گیری زمان اجرای الگوریتم استراسن.

۷. انجام ضرب معمول ماتریس‌ها با استفاده از عملگر داخلی `*` در MATLAB.

۸. اندازه‌گیری زمان اجرای ضرب معمول.

۹. نمایش هر دو نتیجه.

۱۰. محاسبه حداکثر خطای مطلق بین دو نتیجه.

## ۶. آشنایی با پارامتر Cutoff

در این پیاده‌سازی از پارامتری به نام `cutoff` استفاده می‌شود.

برای مثال:

```matlab
n = 10;
cutoff = 2;
```

الگوریتم استراسن به‌صورت بازگشتی ماتریس‌ها را به زیرماتریس‌های کوچک‌تر تقسیم می‌کند.

اما ادامه دادن این تقسیم‌بندی برای ماتریس‌های بسیار کوچک بهینه نیست؛ زیرا فراخوانی‌های بازگشتی و عملیات ماتریسی سربار محاسباتی ایجاد می‌کنند.

پارامتر cutoff مشخص می‌کند که الگوریتم در چه اندازه‌ای تقسیم بازگشتی را متوقف کرده و از ضرب معمول ماتریس‌ها استفاده کند.

برای مثال:

```matlab
cutoff = 2;
```

زمانی که اندازه ماتریس در مرحله بازگشتی به آستانه cutoff برسد، ضرب معمول به‌عنوان حالت پایه (Base Case) اجرا می‌شود.

به این روش، رویکرد ترکیبی (Hybrid Approach) گفته می‌شود؛ زیرا الگوریتم استراسن را با ضرب معمول ماتریس‌ها ترکیب می‌کند.

تغییر مقدار cutoff می‌تواند به ما کمک کند تا بهترین تنظیمات را برای اندازه ماتریس و سیستم موردنظر پیدا کنیم.

## ۷. نحوه اجرای پروژه

### پیش‌نیازها

- MATLAB یا GNU Octave.
- قرار داشتن هر دو فایل `.m` در یک پوشه.

برای اجرای پروژه به هیچ کتابخانه خارجی نیاز نیست.

### اجرای پروژه در GNU Octave

ابتدا ترمینال را باز کنید و وارد پوشه پروژه شوید:

```bash
cd ~/Documents/UNI-semester-3/courseworks/Algo-design-strassen-algo
```

سپس Octave را اجرا کنید:

```bash
octave
```

اسکریپت اصلی را اجرا کنید:

```matlab
strassen_demo
```

همچنین می‌توانید مستقیماً از ترمینال لینوکس اجرا کنید:

```bash
octave --no-gui --quiet strassen_demo.m
```

### اجرای پروژه در MATLAB

MATLAB را باز کنید، وارد پوشه پروژه شوید و دستور زیر را اجرا کنید:

```matlab
strassen_demo
```

**نکته مهم:** هر دو فایل باید در مسیر قابل دسترسی MATLAB یا Octave باشند. نام فایل تابع باید دقیقاً `strassen_multiply.m` باشد؛ زیرا اسکریپت اصلی تابع را با همین نام فراخوانی می‌کند.

## ۸. مقایسه عملکرد دو الگوریتم

در این پروژه دو روش با یکدیگر مقایسه می‌شوند:

**روش اول: الگوریتم استراسن**

```matlab
tic;
C_strassen = strassen_multiply(A, B, cutoff);
t_strassen = toc;
```

**روش دوم: ضرب ماتریس با عملگر داخلی**

```matlab
tic;
C_builtin = A * B;
t_builtin = toc;
```

دستورات `tic` و `toc` برای اندازه‌گیری زمان اجرا استفاده می‌شوند.

برنامه موارد زیر را گزارش می‌کند:

- زمان اجرای الگوریتم استراسن.
- زمان اجرای ضرب داخلی ماتریس‌ها.
- حداکثر خطای مطلق بین دو نتیجه.

### چرا ممکن است ضرب داخلی MATLAB سریع‌تر باشد؟

برای ماتریس‌های کوچک، ضرب داخلی MATLAB ممکن است سریع‌تر باشد؛ زیرا از کتابخانه‌های محاسباتی بهینه‌شده استفاده می‌کند و سربار فراخوانی‌های بازگشتی را ندارد.

مزیت تئوری استراسن بیشتر در اندازه‌های به‌اندازه کافی بزرگ ماتریس‌ها اهمیت پیدا می‌کند.

بنابراین، آزمایش با ماتریس‌های کوچک لزوماً برتری عملکرد استراسن را نشان نمی‌دهد.

برای انجام آزمایش‌های معنادارتر می‌توانید اندازه‌های مختلفی را امتحان کنید:

- 10 × 10
- 100 × 100
- 500 × 500
- 1000 × 1000

همچنین مقدار cutoff را تغییر دهید.

برای مقایسه قابل اعتمادتر زمان اجرا، هر آزمایش را چندین بار تکرار کرده و میانگین زمان‌ها را مقایسه کنید.

## ۹. بررسی دقت محاسباتی (Numerical Accuracy)

در این برنامه حداکثر خطای مطلق با دستور زیر محاسبه می‌شود:

```matlab
max(abs(C_strassen(:) - C_builtin(:)))
```

عملگر `(:)` تمام عناصر ماتریس را به یک بردار ستونی تبدیل می‌کند.

تابع `abs()` اختلاف مطلق عناصر متناظر دو ماتریس را محاسبه می‌کند.

تابع `max()` بزرگ‌ترین اختلاف را پیدا می‌کند.

از نظر ریاضی:

\[
E*{\max} =
\max*{i,j}
|C*{\text{Strassen}}(i,j)-C*{\text{Builtin}}(i,j)|
\]

اگر حداکثر خطا بسیار نزدیک به صفر باشد، دو پیاده‌سازی از نظر عددی نتایج معادل تولید کرده‌اند.

برای مثال:

```text
Maximum Error: 2.22e-15
```

این مقدار اختلاف بسیار کوچکی است که معمولاً در محاسبات عددی ممیز شناور (Floating-Point Arithmetic) انتظار می‌رود.

وجود چنین اختلافی لزوماً به معنی اشتباه بودن الگوریتم نیست.

## ۱۰. خروجی مورد انتظار

اسکریپت اصلی اطلاعات زیر را نمایش می‌دهد:

```text
========== FINAL RESULT ==========

Strassen Result:
[Matrix C]

Built-in Result:
[Matrix C]

Strassen Time: ...
Built-in Time: ...
Maximum Error: ...
```

مقادیر واقعی زمان اجرا و خطای عددی به سیستم کامپیوتری، محیط نرم‌افزاری، اندازه ماتریس و مقدار cutoff بستگی دارند.

## ۱۱. جمع‌بندی و نکات آموزشی

پس از مطالعه و اجرای این پروژه، دانشجو باید بتواند موارد زیر را درک و توضیح دهد:

۱. نحوه انجام ضرب معمول ماتریس‌ها.

۲. چگونگی استفاده از روش تقسیم و حل در ضرب ماتریس‌ها.

۳. دلیل استفاده استراسن از هفت ضرب بازگشتی به‌جای هشت ضرب.

۴. نحوه ترکیب هفت حاصل‌ضرب میانی برای ساخت چهار بلوک ماتریس نهایی.

۵. مفهوم بازگشت (Recursion) و حالت پایه (Base Case).

۶. دلیل استفاده از Padding برای ابعاد ماتریسی که توان عدد ۲ نیستند.

۷. دلیل پیچیدگی زمانی \(O(n^{2.807})\) در الگوریتم استراسن.

۸. دلیل اینکه پیچیدگی زمانی تئوری بهتر همیشه به معنی عملکرد سریع‌تر در تمام اندازه‌ها نیست.

۹. نحوه مقایسه تجربی دو الگوریتم.

۱۰. تأثیر دقت محاسباتی ممیز شناور بر نتایج عددی.
