## Collection of climate-related matplotlib hacks written by Jesse Day.
import matplotlib.pyplot as plt

## Handy function that creates a dual axis with month below and day of year above.
## May wind up having to tweak position of title after running function.
def add_double_calendar_axis():

    ax_ticks = plt.gca()

    ax_ticks.set_xlim([0,365])
    ax_months = ax_ticks.twiny()
    ax_ticks.set_xticks([0,31,59,90,120,151,181,212,243,273,304,334,365])
    ax_ticks.set_xticklabels("")

    ## stealth axis!
    ax_months.set_xlim(ax_ticks.get_xlim())
    ax_months.set_xticks([15.5,45,74.5,105,135.5,165.5,196.5,227.5,258,288.5,319,349.5])
    ax_months.set_xticklabels('JFMAMJJASOND')

    ax_days = ax_months.twiny()

    ax_months.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are off

    ax_days.set_xlim(ax_months.get_xlim())
    ax_days.tick_params(axis='x', which='major', pad=0)
    ax_days.set_xticklabels(["",50,100,150,200,250,300,350])
    
    ## annotate day axis with unit
    ax_ticks.annotate('Day #', (.01, 1.04), textcoords='axes fraction', size=12)

    
## the following draws dashed lines on the plot to demarcate time periods of interest.
## if two time periods are contiguous, then we won't draw two separate lines.
## periods should be a list of lists, with each sublist containing a start and end date.

## also optionally numbers the time periods
def mark_time_periods(periods, numbers = 'no'):
    
    ax = plt.gca()
    
    for i, p in enumerate(periods):
            
        [ymin, ymax] = ax.get_ylim()
        plt.plot([p[0],p[0]], [ymin,ymax],'k--')
        plt.plot([p[1],p[1]], [ymin,ymax],'k--')
        props = dict(boxstyle='square,pad=0.1', facecolor='white')

        if numbers == 'yes':
        #if yes, place a text box in upper left in axes coords
            ax.annotate(i+1, ((p[0]+7)/365, .95), textcoords='axes fraction', \
                        transform=ax.transAxes, fontsize=12,verticalalignment='top', bbox=props)