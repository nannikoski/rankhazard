rankhazard NEWS

CHANGES IN rankhazard 1.1.0:

Visible changes:

* Plenty of new error and warning messages added.
* A bug regarding the use of refpoints in certain 
  situations is fixed.
* Updated examples in the rankhazardplot documentation.
* More flexibility for perfecting the look of the plot
  by new arguments.

Added arguments:

* col.CI
* lty.CI
* lwd.CI
* axes
* xtext
* add
* graphsbefore
* args.legend

Changes in argument names:

* confint -> draw.confint
* refline.col -> col.refline
* refline.lty -> lty.refline
* refline.lwd -> lwd.refline

Non-visible changes:

* The hidden functions coxph_CI and cph_CI are now hidden functions
  rankhazard_CI.coxph, rankhazard_CI.cph and rankhazard_CI and some 
  of the code regarding the calculating of the confidence intervals
  has been moved from rankhazardplot.coxph and rankhazardplot.cph
  into rankhazard_CI.coxph and rankhazard_CI.cph.
* The plotting functions are matplot, matlines and matpoints instead
  of plot, lines and points in the previous version.