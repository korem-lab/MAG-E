import seaborn as sns
import os
from matplotlib import pyplot as plt

PALETTE_KEY = [
    'MaxBin2', 'DAS(MaxBin2-METABAT2-CONCOCT)', 
    'DAS(COMEBin-SemiBin2-METABAT2)', 'VAMB', 'CONCOCT', 'METABAT2', 'SemiBin2', 'COMEBin',
    'multi-sample-all', 'multi-sample-mash', 'single-sample', 'All', 'Prophage', 'Shared'
]
PALETTE = sns.color_palette('Set3', len(PALETTE_KEY))
PALETTE_DICT = dict(zip(PALETTE_KEY, PALETTE))
PALETTE_DICT['multi-sample-all'] = sns.color_palette('Set2', 3)[0]
PALETTE_DICT['multi-sample-mash'] = sns.color_palette('Set2', 3)[1]
PALETTE_DICT['single-sample'] = sns.color_palette('Set2', 3)[2]

MARKER=['o', 'X', 'v', 's', '^']
MARKER_KEY=['multi-sample-all', 'multi-sample-mash', 'single-sample', 'megahit', 'metaspades']
MARKER_DICT=dict(zip(MARKER_KEY, MARKER))

def plot_pipeline_performance_model(emmeans, by, y_name, order):
    
    # filter emmeans
    emmeans = emmeans.loc[emmeans.binner.isin(order),:]
    if by=='assembler':
        strip_hue = 'binner'
        strip_style = 'binning_mode'
        hue_order=order
        order = ['megahit', 'metaspades']
    elif by=='binner':
        strip_hue = 'binning_mode'
        hue_order=['single-sample', 'multi-sample-mash', 'multi-sample-all']
        strip_style = 'assembler'
    elif by=='binning_mode':
        strip_hue = 'binner'
        hue_order = order
        strip_style = 'assembler'
        order=['single-sample', 'multi-sample-mash', 'multi-sample-all']
    else:
        Exception('not valid by')
    ax = None
    for i, cat in enumerate(sorted(emmeans[strip_style].unique())):
        print('marker', cat, MARKER_DICT[cat])
        ax=sns.swarmplot(
            x=emmeans[emmeans[strip_style]==cat][by], 
            y=emmeans[emmeans[strip_style]==cat].emmean, 
            hue=emmeans[emmeans[strip_style]==cat][strip_hue], 
            order=order,
            hue_order=hue_order,
            palette=PALETTE_DICT,
            ax=ax,
            marker=MARKER_DICT[cat],dodge=True,
            s=4.5, linewidth=0.5, edgecolor='black',
            legend=True if i == 0 else False
        )
    ax.set_xlabel('Binners' if by =='binner' else 'Binning Mode' if by =='binning_mode' else 'Assembler')
    ax.set_ylabel(y_name)
    return ax

def plot_pipeline_performance_data(
        score_files, min_rc, min_pr, metric, order, min_cov=2.5
    ):
    fig, ax = plt.subplots(figsize=(8, 2))
    scores, rcvgnms = recoverable_genome_set(score_files, 'gt_recovered', min_cov, min_rc, min_pr)
    scores = scores[scores.metric == metric]
    if metric == 'cov_pr':
        scores.value = 1-scores.value
    scores['pipeline'] = scores.apply(
        lambda x: f'{x.assembler}-{x.binning_mode}', axis=1
    )
    scores = scores.sort_values(by='pipeline')
    ax=sns.boxplot(
        x=scores.pipeline, y=scores.value, showfliers=False, hue=scores.binner,
        ax=ax, palette=PALETTE_DICT, hue_order=order
    )
    return ax

def plot_unlabelled_version(ax, name, tight_layout=True,show=False,ugap=0,lgap=0, ylower=None, yupper=None, xupper=None, xlower=None,remove_top_tick=False, n_ticks=None, xticks=None):
    # first save a fully labeled png version
    if show:
        plt.show()
        return
    if n_ticks:
        from matplotlib.ticker import MaxNLocator
        ax.set_ylim(ylower, yupper)
        ax.xaxis.set_major_locator(MaxNLocator(n_ticks, prune=None))
        ax.yaxis.set_major_locator(MaxNLocator(n_ticks, prune=None))
    else: 
        ymin,ymax = ax.get_ylim()
        xmin,xmax = ax.get_xlim()
        print(ymin-lgap, ymax+ugap)
        ax.set_ylim(ymin-lgap if not ylower else ylower,ymax+ugap if not yupper else yupper)
        ax.set_xlim(xmin if not xlower else xlower, xmax if not xupper else xupper)
        if remove_top_tick:
            yticks = ax.yaxis.get_major_ticks()
            yticks[-1].label1.set_visible(False)
            yticks[-1].tick1line.set_visible(False)
    if xticks:
        ax.set_xticks(xticks)
    plt.savefig(f'{name}_LBLD.png',dpi=900,bbox_inches='tight')
    name = os.path.splitext(name)[0]
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')

    frame1 = plt.gca()
    frame1.legend().set_visible(False)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(f'{name}_UNLBLD.pdf', dpi=1400, bbox_inches='tight')
    plt.savefig(f'{name}_UNLBLD.png', dpi=900, bbox_inches='tight')
    plt.clf()