# LogMoM

This code accompanies the paper "Modeling of bubble coalescence and break-up
using the Log-normal Method of Moments" by V. Habiyaremye, E.M.J. Komen, J.G.M.
Kuerten and E.M.A. Frederix. For a description on what LogMoM does, please refer
to that paper.

## Authors

LogMoM is developed by Edo Frederix and Victor Habiyaremye as part of the
PIONEER programme at the Nuclear Research and Consultancy Group (NRG),
Westerduinweg 3, 1755 LE Petten, the Netherlands. PIONEER is financed by the
Dutch Ministry of Economic Affairs and Climate Policy.

## License

LogMoM is published under the GNU GPL Version 3 license.

LogMoM is distributed under the European Dual Use Codification N: EU DuC=N.
Goods labeled with an EU DuC (European Dual-Use Codification) not equal to 'N'
are subject to European and national export authorization when exported from the
EU and may be subject to national export authorization when exported to another
EU country as well. Even without an EU DuC, or with EU DuC 'N', authorization
may be required due to the final destination and purpose for which the goods are
to be used. No rights may be derived from the specified EU DuC or absence of an
EU DuC.

## Prerequisites

* OpenFOAM-13 Foundation version. While it may compile against other versions,
  this is not tested and currently not supported.
* Python with Numpy, Scipy and Matplotlib for post-processing scripts

## Usage

* Make sure that OpenFOAM-13 is loaded into your environment
* Compile all libraries and apps with

<pre>
./Allwmake
</pre>

* All cases can be run using

<pre>
./prep.sh
foamRun
</pre>

* The prep.sh shell script sets up the case by configuring the required
  parameters. Usually, the prep.sh script takes arguments. If none are
  specified, default arguments will be selected.
* As an example, the TOPFLOW case can also be run in parallel, with

<pre>
./prep.sh
mpirun -np 8 foamRun -parallel
</pre>

* This runs the TOPFLOW case (with default parameters, see prep.sh) in parallel
* Some cases contain a Python post-processing script to plot results. This can
  typically be done with:

<pre>
python plot.py
</pre>

## Using LogMoM in your own OpenFOAM cases

The LogMoM model comes as an additional `diameterModel` which can be selected
for each phase separately. LogMoM is implemented analogously to the standard
IATE diameter model in OpenFOAM. To use LogMoM in your own cases, you must first
compile LogMoM as described above. This will put a shared object called
`libLogMoM.so` in your `$FOAM_USER_LIBBIN` path. Your case must load this shared
object by adding

<pre>
libs
(
    "libLogMoM.so"
);
</pre>

to your case's `system/controlDict` file. Next, you can select the `LogMoM`
model as a `diameterModel` for each phase, in the `constant/phaseProperties`
file. Additional parameters, such as the coalescence and break-up models to use,
must be provided as well. An example is given in
`cases/TOPFLOW/constant/phaseProperties.LogMoM`. Finally, your case should
provide initial and boundary conditions for the `lambda.<phase>` and
`kappai.<phase>` fields, which are the phase `<phase>` zeroth and second
diameter-based size distribution moments scaled by the volume fraction. If a
phase change model is enabled, the `lambda.<phase>` and `kappai.<phase>` fields
must be equipped with the `interfacialGrowthMoment` source, for example, as is
done in the `TOPFLOWCondensation` case. Refer to the cases in `cases/` for all
examples.

## References

* Frederix, E.M.A., Cox, T.L.W., Kuerten, J.G.M., & Komen, E.M.J. (2019).
  Poly-dispersed modeling of bubbly flow using the log-normal size distribution.
  *Chemical Engineering Science*, 201, 237-246.
* Habiyaremye, V., Komen, E.M.J., Kuerten, J.G.M., & Frederix, E.M.A. (2022).
  Modeling of bubble coalescence and break-up using the Log-normal Method of
  Moments. *Chemical Engineering Science*, 253, 117577.
* Habiyaremye, V., Kuerten, J.G.M., & Frederix, E.M.A. (2023). Comparison of
  population balance models for polydisperse bubbly flow. *Chemical Engineering
  Science*, 278, 118932.
* Frederix, E. M. A., Habiyaremye, V., Tajfirooz, S., & Kuerten, J. G. M.
  (2024). Extension of the two-fluid model to bubble size distribution moment
  velocities. *Experimental and Computational Multiphase Flow*, 1-16.

## Contact & support

For bug reports or support, feel free to contact Edo Frederix at
frederix@nrg.eu. Please note that this code is not maintained nor regularly
updated, and is only tested with OpenFOAM-13. Questions related to other
versions will thus not be answered.

## Disclaimer

LogMoM is provided by the copyright holders and contributors "as-is" and any
express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are
disclaimed. In no event shall the copyright owner or contributors be liable for
any direct, indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused and on
any theory of liability, whether in contract, strict liability, or tort
(including negligence or otherwise) arising in any way out of the use of this
software, even if advised of the possibility of such damage.
