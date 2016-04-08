Power Measurement and Savings

Here you will find several addendums to the Power Measurement and
Savings section.

Sample CView files: The following commands will allow you to navigate
through sample data captures to get a feel for what you can do with
CView.

View runs best on Linux; the package can be downloaded from
https://github.com/EMSL-MSC/cview/

Once installed, you can look through the following datasets to
navigate the data; each of these data captures is progressively
larger, with the second and third having more data showing Intel® Xeon
Phi™ coprocessor usage.

cviewall -url http://lotsofcores.com/power/cview1/ -metrics '(aggregation-cpu-average/cpu-user,mic-0/cpu-user,mic-1/cpu-user)’
cviewall -url http://lotsofcores.com/power/cview2/ -metrics '(aggregation-cpu-average/cpu-user,mic-0/cpu-user,mic-1/cpu-user)’

cviewall -url http://lotsofcores.com/power/cview3/ -metrics '(aggregation-cpu-average/cpu-user,mic-0/cpu-user,mic-1/cpu-user)’

CView data navigation instructions can be found here:
https://github.com/EMSL-MSC/cview/wiki/CView-Interaction.  The
tweakbar menus that appear when you start CView will enable you to
increase the size of the data being viewed,

