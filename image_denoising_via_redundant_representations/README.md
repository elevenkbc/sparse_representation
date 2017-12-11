## Image Denoising Via Sparse and Redundant Representation Over Learned Dictionaries 

## 論文重點整理

考慮一個含有噪音的的模型:
$$
y =  x + v   
$$
其中 $y$ 為我們觀察到的訊號 ，$x$ 為真實訊號(未知)，$v$ 為噪音，一般來說我們假設噪音滿足常態分配
$$
v \sim N(0, \sigma)
$$
更進一步假設，我們想要估計的訊號 $x$ ，在某一個字典(Dictionary)  $D$ 中有sparse 的性質，也就是在某一個字典中，$x$ 的 表示係數 $\alpha$ 大多數都是零，也就是說:
$$
x \approx D \alpha
$$
我們稱 $\alpha$  為在字典 $D$ 之下 $x$ 的表示係數。

給定字典 $D$ 與 訊號$x$ 的條件下，求稀疏表示係數 $\alpha$  可以寫成以下最佳化問題:
$$
\left\{ \begin{array}{l} 
& \displaystyle \min_{\alpha} \Vert \alpha\Vert_0\\
& \mbox{subject to} \quad x\approx D\alpha \\
\end{array} \right. \tag{1}
$$
用更確切的符號可以將 $(1)$ 改寫成:
$$
\left\{ \begin{array}{l} 
& \displaystyle \min_{\alpha} \Vert \alpha\Vert_0\\
& \mbox{subject to } \Vert x - D \alpha \Vert_2^2 \leq \epsilon \\
\end{array} \right. \tag{2}
$$
求稀疏表示係數 $\alpha$  也可以寫成以下最佳化問題:
$$
\left\{ \begin{array}{l} 
& \displaystyle \min_{\alpha} \Vert x - D \alpha \Vert_2^2\\
& \mbox{subject to } \Vert  \alpha \Vert_0 \leq L \quad L \ll \dim(x)\\
\end{array} \right. \tag{3}
$$
其中 $\Vert \alpha \Vert_0$ 為 L0範數，表示 $\alpha$向量中的非零元的個數，$\epsilon$ 為 $x$ 訊號的稀疏表示誤差。

問題$(2)$與問題$(3)$ 非常相似，一般來說$(2)$與$(3)$都是NP hard ，因為L0 範數涉及到組合搜尋問題。

不過有一些不錯的數值方法可以來求$(2)$與 $(3)$的近似解，例如OMP ，

在以上這些演算法中，又屬OMP 最簡單易實現，計算效率也高，故在這裡提出的演算法的sparse coding 步驟都是基於 OMP演算法。

接著我們來討論如何用給定的字典 $D$ 來對影像做去噪，其中字典 $D$ 能將影像的每個分塊影像作稀疏表示。

### 分塊影像的稀疏模型(sparseland Mode)

#### 定義 $(\epsilon, L, D) - \mbox{Sparseland signals}$:

若 $x$ 訊號對於給定 Triple $(\epsilon, L, D)$ 存在一個向量 $\alpha$滿足:
$$
\Vert  D \alpha - x\Vert_2^2 \leq \epsilon, \quad \Vert \alpha \Vert_0 \leq L \ll \dim(x)
$$
稱 $x$ 為 $(\epsilon, L, D) - \mbox{Sparseland signals}$ . 

 影像 $I \in \mathbf{R}^{\sqrt{N} \times \sqrt{N}}$ ，假設成$\sqrt{N}$ 的理由是讓 $I$ 可以向量化成 $\mathbf{R}^{N \times 1}$ 的行向量。

假設受雜訊污染的影像 $I_{noise}$由以下方程式給出:
$$
I_{noise} = I +  \mbox{noise} \tag{4}
$$
我們假設雜訊為白噪音，更確切的說式 $(4)$ 可寫成:
$$
(I_{noise})_{ij} = (I)_{ij} + v  \quad  \mbox{ for i, j = } 1,2, ...,\sqrt{N}
$$
其中 $v \sim N(0, \sigma)$.

為了方必定義提取分塊影像的算子，我們將影像都拉成行向量來考慮，令 $X  := \mbox{Vec} (I)$  和 $Y := \mbox{Vec}(I_{noise})$ 和 $V := \mbox{Vec}(\mbox{noise})$ ， $X , Y, V \in \mathbf{R}^{N \times 1}$ 。

其中 $\mbox{Vec}()$ 為矩陣向量化運算，定義如下:

$\mbox{Vec}(A) = [A_1^T, A_2^T, ... A_N^T]^T$ ，其中 $A_i$ 為 $A$矩陣中的第 $i$行。

因此 $(4)$ 式可以改寫如下:
$$
Y = X + V
$$


將拉成行向量的影像 $X$ 分成許多小塊的分塊影像 $x_{ij} \in \mathbf{R^{n\times 1}}$ ，這些分塊影像 $x_{ij}$ 在原始影像 $I$ 中可能是互相重疊的區塊。

定義提取分塊影像矩陣 $R_{ij}  \in \mathbf{R}^{n \times N}$，使得 $R_{ij}X = x_{ij}$，其中 $i, j $ 分別代表分塊影像矩陣的列索引與行索引。

假設每個分塊影像都滿足 $x_{ij}$ 均是 $(\epsilon, L, D) - \mbox{Sparseland signals}$ ，則對於每一個分塊影像 $x_{ij}$，可以利用以下模型來描述字典 $D$ 的稀疏表示問題:
$$
\hat{\alpha}_{ij} = \arg \min_\alpha \Vert \alpha \Vert_0  \quad  \mbox{   subject to   }  \Vert D \alpha  - x_{ij}\Vert_2^2 \leq \epsilon\tag{5}
$$

但是在影像去噪的問題中，我們並沒有真實的分塊影像 $x_{ij}$ ，我們只有受雜訊干擾的分塊影像 $y_{ij}$ ，$y_{ij}$定義如下:
$$
y_{ij} = R_{ij}Y =  R_{ij}(X +V) = R_{ij}X + R_{ij}V= x_{ij} + v_{ij}
$$
注意 $y_{ij}, x_{ij}, v_{ij} \in \mathbf{R}^{n \times 1}$ ，在未知真實的分塊影像 $x_{ij} $ 的前提下，問題$(5)$需改寫為:
$$
\hat{\alpha}_{ij} = \arg \min_\alpha \Vert \alpha \Vert_0  \quad  \mbox{   subject to   }  \Vert D \alpha - y_{ij} \Vert \leq T \tag{6}
$$
這裡的 $T$ 與雜訊標準差 $\sigma$ 和 稀疏表示誤差 $\epsilon$ 有關。

注意，由最佳化理論我們知道，存在一個  $\mu > 0$ 使得問題$(6)$ 等價於以下無限制最佳化問題:
$$
\hat{\alpha}_{ij} = \arg \min_\alpha \Vert D\alpha - y_{ij} \Vert_2^2 + \mu \Vert \alpha \Vert_0   \tag{7}
$$
對於每一個受雜訊干擾的分塊影像 $y_{ij}$，我們可以透過利用 OMP 求解 問題$(7)$  ，接著利用  $\hat{\alpha}_{ij}$ 與字典 $D$ 來估計分塊影像 $x_{ij}$，估計的公式由以下方程式給出:
$$
\hat{x}_{ij} = D \hat{\alpha}_{ij}
$$
下一個節，我們討潤大影像中的去躁問題，並寫出大影像去噪問題相對應的最佳化模型。

## 分塊影像的去噪音到大影像的去噪

假設大影像 $X$ 中的每個分塊影像都滿足 $x_{ij}$ 均是 $(\epsilon, L, D) - \mbox{Sparseland signals}$，且受雜訊干擾的圖片 $Y$ 與原圖 $X$ 滿足以下關係:
$$
\Vert X - Y \Vert _2^2 \leq \mbox{Constant} \cdot \sigma^2 \tag{8}
$$
其中 $\mbox{Constant}$ 與影像維度 $N$ 有關，且 $\sigma$ 為白噪音 $v$ 的標準差 。

$\Vert X - Y \Vert _2^2 $ 是一個有界的項，因此我們可以將式 $(8)$ 當作一個逞罰項加入 問題 $(5)$ 中，由最佳化理論可以將大影像的去躁問題可以改寫成以下無限制最佳化問題:
$$
\{ \hat{\alpha}_{ij}, \hat{X}\} =  \arg \min_{\alpha_{ij}, X} \lambda \Vert Y-X\Vert_2^2+ \sum_{ij}\mu_{ij} \Vert \alpha_{ij}\Vert_0 + \sum_{ij} \Vert D\alpha_{ij} - R_{ij}X\Vert_2^2 \tag{9}
$$

其中  $\lambda $為與 $\sigma$ 有關的參數，且 $\lambda > 0$。

## 數值解法

最佳化問題 $(9)$ 中，假設稀疏表示字典 $D$ 已知，問題 $(9)$ 中想要求解的變數有兩個，一個是去造後的影像 $\hat{X}$ ，一個是每個分塊矩陣的稀疏表示係數 $\hat{\alpha}_{ij}$ ，同時求解這兩個變數是很困難的，故本論文提出了兩階段的方法來求解問題 $(9)$，首先假設 $Y = X$ 故 $(9)$ 可以化簡為:
$$
\hat{\alpha}_{ij} =  \arg \min_{\alpha_{ij}, X}  \sum_{ij}\mu_{ij} \Vert \alpha_{ij}\Vert_0 + \sum_{ij} \Vert D\alpha_{ij} - R_{ij}X\Vert_2^2 \tag{10}
$$
我們可以透過 OMP 算法來求解 $(9)$，得到每個分塊矩陣的稀疏表示係數 $\hat{\alpha}_{ij}$，在已知的  $\hat{\alpha}_{ij}$ 的情況下問題 $(9)$ 可以化簡為:
$$
\hat{X} =  \arg \min_{\alpha_{ij}, X} \lambda \Vert Y-X\Vert_2^2 + \sum_{ij} \Vert D\alpha_{ij} - R_{ij}X\Vert_2^2 \tag{11}
$$
最佳化問題 $(10)$ 有close-form解，令
$$
\begin{align*}
L(X) &= \lambda \Vert Y-X\Vert_2^2 + \sum_{ij} \Vert D\alpha_{ij} - R_{ij}X\Vert_2^2\\
& =  \lambda ( Y-X)^T( Y-X) +  \sum_{ij}  (D\alpha_{ij} - R_{ij}X)^T(D\alpha_{ij} - R_{ij}X) \\
& =\lambda Y^TY - 2 \lambda X^TY +\lambda X^TX +  \sum_{ij} \alpha_{ij}^TD^TD\alpha_{ij} - 2X^TR_{ij}^TD\alpha_{ij} + X^TR_{ij}^TR_{ij}X
\end{align*}
$$
$L(X)​$ 對向量 $X​$ 微分:
$$
\begin{align*}
\dfrac{\partial }{\partial X}L(X) &= - 2\lambda Y + 2\lambda X +  \sum_{ij}  -2R_{ij}^TD\alpha_{ij}  + 2 R_{ij}^TR_{ij}X
\end{align*}
$$
令 $  \dfrac{\partial }{\partial X}L(X)  = 0$
$$
\begin{align*}
 & -2\lambda Y + 2\lambda X +  \sum_{ij}  -2R_{ij}^TD\alpha_{ij}  + 2 R_{ij}^TR_{ij}X = 0 \\
 & -\lambda Y + \lambda X +  \sum_{ij}  -R_{ij}^TD\alpha_{ij}  +  R_{ij}^TR_{ij}X = 0 \\
 &  \lambda X +  \sum_{ij}   R_{ij}^TR_{ij}X = \lambda Y + \sum_{ij}  R_{ij}^TD\alpha_{ij} \\
 & (\lambda I + \sum_{ij}   R_{ij}^TR_{ij})X =  \lambda Y + \sum_{ij}  R_{ij}^TD\alpha_{ij} \\
&  X =  (\lambda I + \sum_{ij}   R_{ij}^TR_{ij})^{-1}( \lambda Y + \sum_{ij}  R_{ij}^TD\alpha_{ij} )
\end{align*}
$$
故我們得到 問題$(10)$的 close-form 為
$$
\hat{X} =  (\lambda I+ \sum_{ij}   R_{ij}^TR_{ij})^{-1}( \lambda Y + \sum_{ij}  R_{ij}^TD\alpha_{ij} )
$$


