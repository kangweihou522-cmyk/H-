"""Draw original teaching diagrams with synthetic numbers, never analysis results."""
from pathlib import Path
import os
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle

OUT = ROOT/'docs/assets'
GREEN, INK, BLUE, RED, GRAY = '#126657', '#183531', '#4275ad', '#d76b59', '#87958f'


def setup():
    candidates = ['Noto Sans CJK SC', 'Microsoft YaHei', 'SimHei', 'Noto Sans CJK JP']
    available = {font.name for font in font_manager.fontManager.ttflist}
    selected = next((name for name in candidates if name in available), None)
    if not selected:
        raise SystemExit('生成教学图需要中文字体；已提交的图片无需在用户电脑重新生成。')
    plt.rcParams.update({'font.family': selected, 'font.size': 16, 'axes.unicode_minus': False,
                         'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': GRAY,
                         'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
    OUT.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.savefig(OUT/name, dpi=160, bbox_inches='tight', metadata={'Software': 'RNA-seq workstation educational diagrams'})
    plt.close(fig)


def arrow(ax, start, end, color=GREEN):
    ax.annotate('', xy=end, xytext=start, arrowprops={'arrowstyle': '-|>', 'color': color, 'lw': 2.2})


def box(ax, x, y, w, h, text, color='#e7f2ec', fontsize=16):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=.012', facecolor=color, edgecolor='#b8cfc2'))
    ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=fontsize)


def molecules():
    fig, ax = plt.subplots(figsize=(12, 7.2))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('一个基因可以产生多个转录本；双端读段来自同一个文库片段', fontsize=21, pad=20)
    colors = [GREEN, BLUE, RED]
    ax.text(.02, .84, '基因组DNA', va='center')
    ax.plot([.25, .93], [.84, .84], color=GRAY, lw=2)
    for x, label, color in zip([.28, .51, .76], ['外显子1', '外显子2', '外显子3'], colors):
        ax.add_patch(Rectangle((x, .805), .14, .07, color=color))
        ax.text(x+.07, .84, label, color='white', fontsize=14, ha='center', va='center')
    ax.text(.46, .895, '内含子', fontsize=12, ha='center')
    ax.text(.71, .895, '内含子', fontsize=12, ha='center')
    arrow(ax, (.6, .78), (.6, .70))
    ax.text(.66, .73, '转录、剪接', fontsize=14)
    for y, label, exons in [(.60, '转录本A', [0, 1, 2]), (.45, '转录本B', [0, 2])]:
        ax.text(.02, y+.035, label, va='center')
        for index, exon in enumerate(exons):
            x = .32+index*.14
            ax.add_patch(Rectangle((x, y), .14, .07, color=colors[exon]))
            ax.text(x+.07, y+.035, str(exon+1), color='white', ha='center', va='center')
    ax.text(.79, .49, '同一基因\n不同剪接形式', fontsize=14, va='center')
    arrow(ax, (.51, .42), (.51, .34))
    box(ax, .03, .12, .31, .15, '片段化、反转录、加接头\n形成cDNA文库', fontsize=15)
    arrow(ax, (.36, .195), (.44, .195))
    ax.plot([.47, .95], [.20, .20], color=GRAY, lw=10, solid_capstyle='round')
    arrow(ax, (.47, .25), (.65, .25), BLUE)
    arrow(ax, (.95, .15), (.77, .15), RED)
    ax.text(.56, .30, 'R1 →', color=BLUE, ha='center', fontsize=16)
    ax.text(.86, .08, '← R2', color=RED, ha='center', fontsize=16)
    ax.text(.71, .36, '同一个双端文库片段', ha='center', fontsize=15)
    ax.text(.5, -.01, '示意图，不按真实长度比例绘制；一对R1/R2不是两个生物学重复。', ha='center', fontsize=13)
    save(fig, 'rna-to-reads.png')


def normalization():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    for ax, values, title, ylabel in zip(axes, [[100, 200], [100, 100]],
            ['原始计数：测得更多不一定表达更高', '按文库大小缩放的教学示例'], ['该基因的计数', '示例CPM（每百万计数）']):
        bars = ax.bar(['样本C\n总量100万', '样本T\n总量200万'], values, color=[BLUE, RED], width=.55)
        ax.set(title=title, ylabel=ylabel, ylim=(0, 250))
        ax.spines[['top', 'right']].set_visible(False)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x()+bar.get_width()/2, value+7, str(value), ha='center')
    fig.suptitle('若表达组成相同，测序深度翻倍可以让计数翻倍', fontsize=22)
    fig.text(.5, .015, '100 / 1,000,000 × 1,000,000 = 100；200 / 2,000,000 × 1,000,000 = 100\n此处为CPM算术示意；DESeq2使用更稳健的size factor及相应长度信息，并非直接做CPM检验。', ha='center', fontsize=12)
    fig.tight_layout(rect=(0, .13, 1, .92))
    save(fig, 'library-normalization.png')


def thresholds():
    genes = np.array(['A', 'B', 'C', 'D', 'E', 'F', 'H'])
    lfc = np.array([.4, 1.4, -1.8, 2.5, 1., 0., -.8])
    p = np.array([.04, .001, .01, .2, .003, 1., .01])
    padj = np.array([.08, .01, .04, .3, .02, 1., np.nan])
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4))
    for ax, column, cutoff, title, label in zip(axes, [p, padj], [0, 1],
            ['方案一：原始P＜0.05，|LFC|＞0', '方案二：padj＜0.05，|LFC|＞1'], ['原始P值', '校正P值padj']):
        selected = (column < .05) & (np.abs(lfc) > cutoff)
        colors = [RED if ok and x > 0 else BLUE if ok else GRAY for ok, x in zip(selected, lfc)]
        ax.scatter(lfc, -np.log10(column), c=colors, s=100, edgecolors=INK, linewidths=.5)
        ax.axhline(-np.log10(.05), color=GRAY, ls='--')
        for x in sorted(set([-cutoff, cutoff])):
            ax.axvline(x, color=GRAY, ls='--')
        for gene, x, y in zip(genes, lfc, column):
            if np.isfinite(y): ax.annotate(gene, (x, -np.log10(y)), xytext=(6, 8), textcoords='offset points')
        ax.set(title=title, xlabel='log2（处理 / 对照）', ylabel='-log10（'+label+'）', xlim=(-2.5, 3), ylim=(-.15, 3.5))
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.13)
    fig.suptitle('同一张结果表，不同筛选规则得到不同候选列表', fontsize=22)
    fig.text(.5, .025, '纯教学示例，不是实验结果。红：入选上调；蓝：入选下调；灰：未入选。\nE的LFC恰为1，不满足严格“＞1”；H的padj缺失，因此不在右图强行补点。', ha='center', fontsize=12)
    fig.tight_layout(rect=(0, .12, 1, .92))
    save(fig, 'deg-two-profiles.png')


def methods():
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.set_title('ORA、GSEA、ssGSEA：三个不同的问题', fontsize=23, pad=15)
    rows = [(.71, 'ORA：看选中的名单', '输入：差异基因列表 + 合适背景\n问题：某功能的成员是否比预期更多？\n例：预期命中4个，实际命中12个', '#e7f2ec'),
            (.40, 'GSEA：看全体基因排序', '输入：带正负方向的全基因排序\n问题：某集合的成员是否集中在排序两端？\n输出：ES、NES、校正P、leading-edge', '#edf2f8'),
            (.09, 'ssGSEA：逐个样本打分', '输入：基因 × 样本表达矩阵 + 基因集\n问题：该样本中，集合成员是否相对靠前？\n先得分数，再用模型比较组别', '#fcf0ea')]
    for y, title, text, color in rows:
        box(ax, .03, y, .30, .22, title, color, 18)
        arrow(ax, (.35, y+.11), (.43, y+.11))
        box(ax, .45, y, .52, .22, text, color, 16)
    ax.text(.5, -.01, '通路富集、转录层面评分都不等于已测得蛋白活性、代谢通量或因果机制。', ha='center', fontsize=13)
    save(fig, 'enrichment-methods.png')


if __name__ == '__main__':
    setup()
    for draw in [molecules, normalization, thresholds, methods]:
        draw()
    print('Generated four original educational diagrams with synthetic data.')
