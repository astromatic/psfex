<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
        <!ENTITY Beta "&#946;">
        <!ENTITY chi "&#967;">
        <!ENTITY darr "&#8595;">
	]>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<!-- 
#				psfex.xsl
#
# Global XSL template
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	PSFEx
#
#	Copyright:		(C) 2005-2012 Chiara Marmo - IAP/CNRS
#				& Emmanuel Bertin - IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	PSFEx is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	PSFEx is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with PSFEx. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		17/06/2012
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -->

 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <HTML>
   <HEAD>
    <link rel="shortcut icon" type="image/x-icon" href="http://astromatic.net/xsl/favicon.ico" />
    <script type="text/javascript" src="http://astromatic.net/xsl/sorttable.js"/>

    <style type="text/css">
     p {
      font-family: sans-serif;
      }
     p.italic {font-style: italic}
     body {
      margin: 10px;
      background-color: #e0e0e0;
      background-image: url("http://astromatic.net/xsl/body_bg.jpg");
      background-repeat: repeat-x;
      background-position: top;
      min-width:662px;
      }
     mono {font-family: monospace}
     elmin {color: green }
     elmax {color: red }
     el {
      font-family: monospace;
      font-size: 100%;
      color: black;
      }
     elm {
      font-family: monospace;
      font-size: 67%;
      white-space: nowrap;
      }
     a {text-decoration: none; font-style: bold; color: #476674}
     a:hover {text-decoration: underline;}
     #header {
      padding: 5px;
      min-width: 662px;
      background-image: url("http://astromatic.net/xsl/astromaticleft.png");
      background-repeat: repeat-x;
      background-position: left top;
      text-align: left;
      font-size: 1.2em;
      margin: 0 0 30px 0;
      color:#d3e7f0;
      font-weight: bold;
      }
     th {
      background-color:#d3e7f0;
      border-top: 1px solid white;
      border-left: 1px solid white;
      border-right: 1px solid #476674;
      border-bottom: 1px solid #476674;
      -moz-border-radius: 3px;
      -khtml-border-radius: 3px;
      -webkit-border-radius: 3px;
      border-radius: 3px;
      padding: 2px;
      line-height: 12px;
      }
     td {
      background-color:#f2f4f4;
      padding-left: 2px;
      padding-right: 2px;
      }
     table.sortable {
      border-top: 1px solid #476674;
      border-left: 1px solid #476674;
      border-right: 1px solid white;
      border-bottom: 1px solid white;
      -moz-border-radius: 3px;
      -khtml-border-radius: 3px;
      -webkit-border-radius: 3px;
      border-radius: 3px;
      }
     table.sortable a.sortheader {
      background-color:#d3e7f0;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }

     table.sortable span.sortarrow {
      color: black;
      font-weight: bold;
      text-decoration: blink;
      }
     table.sortable a.sortheader.sub {vertical-align: sub}
     </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </HEAD>
    <BODY>
     <div id="header">
      <a href="http://astromatic.net"><img style="vertical-align: middle; border:0px" src="http://astromatic.net/xsl/astromatic.png" title="Astromatic home" alt="Astromatic.net" /></a>  Processing summary
     </div>
    <xsl:call-template name="VOTable"/>
   </BODY>
  </HTML>
 </xsl:template>

<!-- **************** Generic XSL template for VOTables ****************** -->
 <xsl:template name="VOTable">
  <xsl:for-each select="/VOTABLE">
   <xsl:call-template name="Resource"/>
  </xsl:for-each>
 </xsl:template>

<!-- *************** Generic XSL template for Resources ****************** -->
 <xsl:template name="Resource">
  <xsl:for-each select="RESOURCE">
   <xsl:choose>
    <xsl:when test="@ID='PSFEx'">
     <xsl:call-template name="PSFEx"/>
    </xsl:when>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>

<!-- ********************** XSL template for PSFEx *********************** -->
 <xsl:template name="PSFEx">
  <xsl:for-each select="RESOURCE[@ID='MetaData']">
   <xsl:call-template name="RunInfo"/>
   <xsl:for-each select="TABLE[@ID='PSF_Fields']">
    <xsl:call-template name="PSF_Fields"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='PSF_Extensions']">
    <xsl:call-template name="PSF_Extensions"/>
   </xsl:for-each>
   <xsl:for-each select="RESOURCE[@ID='Config']">
    <xsl:call-template name="Config"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Warnings']">
    <xsl:call-template name="Warnings"/>
   </xsl:for-each>
  </xsl:for-each>
 </xsl:template>

<!-- ************* Generic XSL RunInfo template for MetaData ************* -->
 <xsl:template name="RunInfo">
  <p>
<!-- Software name, version, date, time and number of threads -->
   <a>
    <xsl:attribute name="href">
     <xsl:value-of select="PARAM[@name='Soft_URL']/@value"/>
    </xsl:attribute>
    <b>
     <xsl:value-of select="PARAM[@name='Software']/@value"/>&nbsp;<xsl:value-of select="PARAM[@name='Version']/@value"/>
    </b>
   </a>
   started on
   <b><xsl:value-of select="PARAM[@name='Date']/@value"/></b>
   at
   <b><xsl:value-of select="PARAM[@name='Time']/@value"/></b>
   with
   <b><xsl:value-of select="PARAM[@name='NThreads']/@value"/></b>
   thread<xsl:if test="PARAM[@name='NThreads']/@value &gt; 1">s</xsl:if>

<!-- Run time -->
   <xsl:variable name="duration" select="PARAM[@name='Duration']/@value"/>
   (run time:
    <b>
     <xsl:choose> 
      <xsl:when test="$duration &gt; 3600.0">
       <xsl:value-of
	select='concat(string(floor($duration div 3600)),
	" h ", format-number(floor(($duration div 60) mod 60.0), "00"),
	" min")'/>
      </xsl:when>
      <xsl:otherwise>
       <xsl:choose>
        <xsl:when test="$duration &gt; 60.0">
         <xsl:value-of
	  select='concat(format-number(floor($duration div 60),"##"),
	  " min ", format-number(floor($duration mod 60.0), "00")," s")'/>
        </xsl:when>
        <xsl:otherwise>
         <xsl:value-of select='concat(string($duration), " s")'/>
        </xsl:otherwise>
       </xsl:choose>
      </xsl:otherwise>
     </xsl:choose>
    </b>)
    <br />
   by user <b><xsl:value-of select="PARAM[@name='User']/@value"/></b>
   from <b><xsl:value-of select="PARAM[@name='Host']/@value"/></b>
   in <b><mono><xsl:value-of select="PARAM[@name='Path']/@value"/></mono></b>
  </p>
  <p>
   <b style="color: red"><xsl:if test="PARAM[@name='Error_Msg']/@value &gt; 0">
    An Error occured!!! </xsl:if>
   <xsl:value-of select="PARAM[@name='Error_Msg']/@value"/></b>
  </p>
 </xsl:template>

<!-- ********************** XSL template for PSF_Fields ******************* -->
  <xsl:template name="PSF_Fields">
   <xsl:variable name="name" select="count(FIELD[@name='Catalog_Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ident" select="count(FIELD[@name='Image_Ident']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="next" select="count(FIELD[@name='NExtensions']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nloadmin" select="count(FIELD[@name='NStars_Loaded_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nload" select="count(FIELD[@name='NStars_Loaded_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nloadmax" select="count(FIELD[@name='NStars_Loaded_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="naccepmin" select="count(FIELD[@name='NStars_Accepted_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="naccep" select="count(FIELD[@name='NStars_Accepted_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="naccepmax" select="count(FIELD[@name='NStars_Accepted_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fluxradmin" select="count(FIELD[@name='FWHM_FromFluxRadius_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fluxrad" select="count(FIELD[@name='FWHM_FromFluxRadius_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fluxradmax" select="count(FIELD[@name='FWHM_FromFluxRadius_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="samplingmin" select="count(FIELD[@name='Sampling_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="sampling" select="count(FIELD[@name='Sampling_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="samplingmax" select="count(FIELD[@name='Sampling_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2min" select="count(FIELD[@name='Chi2_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2" select="count(FIELD[@name='Chi2_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2max" select="count(FIELD[@name='Chi2_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhmmin" select="count(FIELD[@name='FWHM_WCS_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhm" select="count(FIELD[@name='FWHM_WCS_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhmmax" select="count(FIELD[@name='FWHM_WCS_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipticitymin" select="count(FIELD[@name='Ellipticity_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipticity" select="count(FIELD[@name='Ellipticity_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipticitymax" select="count(FIELD[@name='Ellipticity_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="betamin" select="count(FIELD[@name='MoffatBeta_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="beta" select="count(FIELD[@name='MoffatBeta_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="betamax" select="count(FIELD[@name='MoffatBeta_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="residualsmin" select="count(FIELD[@name='Residuals_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="residuals" select="count(FIELD[@name='Residuals_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="residualsmax" select="count(FIELD[@name='Residuals_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pffwhmmin" select="count(FIELD[@name='FWHM_PixelFree_WCS_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pffwhm" select="count(FIELD[@name='FWHM_PixelFree_WCS_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pffwhmmax" select="count(FIELD[@name='FWHM_PixelFree_WCS_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfellipticitymin" select="count(FIELD[@name='Ellipticity_PixelFree_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfellipticity" select="count(FIELD[@name='Ellipticity_PixelFree_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfellipticitymax" select="count(FIELD[@name='Ellipticity_PixelFree_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfbetamin" select="count(FIELD[@name='MoffatBeta_PixelFree_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfbeta" select="count(FIELD[@name='MoffatBeta_PixelFree_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfbetamax" select="count(FIELD[@name='MoffatBeta_PixelFree_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfresidualsmin" select="count(FIELD[@name='Residuals_PixelFree_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfresiduals" select="count(FIELD[@name='Residuals_PixelFree_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfresidualsmax" select="count(FIELD[@name='Residuals_PixelFree_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="asymmetrymin" select="count(FIELD[@name='Asymmetry_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="asymmetry" select="count(FIELD[@name='Asymmetry_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="asymmetrymax" select="count(FIELD[@name='Asymmetry_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="noiseqareamin" select="count(FIELD[@name='Area_Noise_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="noiseqarea" select="count(FIELD[@name='Area_Noise_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="noiseqareamax" select="count(FIELD[@name='Area_Noise_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="countsplotflag" select="count(FIELD[@name='Plot_Counts'])"/>
   <xsl:variable name="countsplot" select="count(FIELD[@name='Plot_Counts']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="countfracplotflag" select="count(FIELD[@name='Plot_Count_Fraction'])"/>
   <xsl:variable name="countfracplot" select="count(FIELD[@name='Plot_Count_Fraction']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhmplotflag" select="count(FIELD[@name='Plot_FWHM'])"/>
   <xsl:variable name="fwhmplot" select="count(FIELD[@name='Plot_FWHM']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipplotflag" select="count(FIELD[@name='Plot_Ellipticity'])"/>
   <xsl:variable name="ellipplot" select="count(FIELD[@name='Plot_Ellipticity']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" onclick="showhideTable('psffields')" title="click to expand">
     PSF stats per Input File&nbsp;&darr;
    </BUTTON>
    <TABLE class="sortable" id="psffields" BORDER="2" style="display: none">
     <TR>
      <TH>Filename</TH>
      <TH>Identifier</TH>
      <TH>N<sub>ext</sub></TH>
      <TH>N<sub>loaded</sub></TH>
      <xsl:if test="$countsplotflag &gt; 0">
       <TH>Count map</TH>
      </xsl:if>
      <TH>N<sub>accepted</sub></TH>
      <xsl:if test="$countfracplotflag &gt; 0">
       <TH>Count fraction map</TH>
      </xsl:if>
      <TH>Half-Light diam.</TH>
      <TH>Sampling</TH>
      <TH>&chi;<sup>2</sup>/d.o.f.</TH>
      <TH>FWHM</TH>
      <TH>FWHM<BR/>
       <elm>pixel free</elm>
      </TH>
      <xsl:if test="$fwhmplotflag &gt; 0">
       <TH>FWHM map</TH>
      </xsl:if>
      <TH>Ellipticity</TH>
      <xsl:if test="$ellipplotflag &gt; 0">
      <TH>Ellipticity<BR/>
       <elm>pixel free</elm>
      </TH>
       <TH>Ellipticity map</TH>
      </xsl:if> 
      <TH style="white-space: nowrap">Moffat &Beta;</TH>
      <TH style="white-space: nowrap">Moffat &Beta;<BR/>
       <elm>pixel free</elm>
      </TH>
      <TH>Residuals</TH>
      <TH>Residuals<BR/>
       <elm>pixel free</elm>
      </TH>
      <TH>Asymmetry</TH>
      <TH>Noise Equivalent Area</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
<!-- Catalog name -->
        <td>
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
<!-- Identifier -->
        <td align="center">
         <el><xsl:value-of select="TD[$ident]"/></el>
        </td>
<!-- Number of extensions -->
        <td align="center">
         <el><xsl:value-of select="TD[$next]"/></el>
        </td>
<!-- Number of sources loaded -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$nload],'######0.0')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="TD[$nloadmin]"/></elmin>
	  - <elmax><xsl:value-of select="TD[$nloadmax]"/></elmax>
         </elm>
        </td>
<!-- count fraction map -->
        <xsl:if test="$countsplotflag &gt; 0">
         <td align="center">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$countsplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$countsplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$countsplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
<!-- Number of sources accepted -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$naccep],'######0.0')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="TD[$naccepmin]"/></elmin>
	  - <elmax><xsl:value-of select="TD[$naccepmax]"/></elmax>
         </elm>
        </td>
<!-- fraction map -->
        <xsl:if test="$countfracplotflag &gt; 0">
         <td align="center">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$countfracplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$countfracplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$countfracplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
<!-- Half-light diameter -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$fluxrad],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$fluxradmin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$fluxradmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Sampling -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$sampling],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$samplingmin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$samplingmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Chi2 --> 
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$chi2], '######0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$chi2min],'######0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$chi2max],'######0.00')"/></elmax>
         </elm>
        </td>
<!-- FWHM -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$fwhm],'##0.00')"/>"</el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$fwhmmin],'##0.00')"/>"</elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$fwhmmax],'##0.00')"/>"</elmax>
         </elm>
        </td>
<!-- FWHM pixel-free-->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pffwhm],'##0.00')"/>"</el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pffwhmmin],'##0.00')"/>"</elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$pffwhmmax],'##0.00')"/>"</elmax>
         </elm>
        </td>
<!-- FWHM map -->
        <xsl:if test="$fwhmplotflag &gt; 0">
         <td align="center">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$fwhmplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$fwhmplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$fwhmplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
<!-- ellipticity -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$ellipticity],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$ellipticitymin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$ellipticitymax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- ellipticity pixel-free -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pfellipticity],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pfellipticitymin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$pfellipticitymax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Ellipticity map -->
        <xsl:if test="$ellipplotflag &gt; 0">
         <td align="center">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$ellipplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$ellipplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$ellipplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
<!-- Beta -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$beta],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$betamin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$betamax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Beta pixel-free -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pfbeta],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pfbetamin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$pfbetamax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Residuals -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$residuals],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$residualsmin],'##0.00')"/></elmin>
 	  - <elmax><xsl:value-of select="format-number(TD[$residualsmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Residuals pixel-free -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pfresiduals],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pfresidualsmin],'##0.00')"/></elmin>
 	  - <elmax><xsl:value-of select="format-number(TD[$pfresidualsmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Asymmetry -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$asymmetry],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$asymmetrymin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$asymmetrymax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Noise Equivalent Area -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$noiseqarea],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$noiseqareamin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$noiseqareamax],'##0.00')"/></elmax>
         </elm>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

<!-- ********************** XSL template for PSF_Extensions ******************* -->
  <xsl:template name="PSF_Extensions">
   <xsl:variable name="ext" select="count(FIELD[@name='Extension']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nloadmin" select="count(FIELD[@name='NStars_Loaded_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nload" select="count(FIELD[@name='NStars_Loaded_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nloadmax" select="count(FIELD[@name='NStars_Loaded_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="naccepmin" select="count(FIELD[@name='NStars_Accepted_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="naccep" select="count(FIELD[@name='NStars_Accepted_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="naccepmax" select="count(FIELD[@name='NStars_Accepted_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fluxradmin" select="count(FIELD[@name='FWHM_FromFluxRadius_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fluxrad" select="count(FIELD[@name='FWHM_FromFluxRadius_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fluxradmax" select="count(FIELD[@name='FWHM_FromFluxRadius_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="samplingmin" select="count(FIELD[@name='Sampling_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="sampling" select="count(FIELD[@name='Sampling_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="samplingmax" select="count(FIELD[@name='Sampling_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2min" select="count(FIELD[@name='Chi2_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2" select="count(FIELD[@name='Chi2_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2max" select="count(FIELD[@name='Chi2_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhmmin" select="count(FIELD[@name='FWHM_WCS_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhm" select="count(FIELD[@name='FWHM_WCS_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhmmax" select="count(FIELD[@name='FWHM_WCS_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipticitymin" select="count(FIELD[@name='Ellipticity_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipticity" select="count(FIELD[@name='Ellipticity_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipticitymax" select="count(FIELD[@name='Ellipticity_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="betamin" select="count(FIELD[@name='MoffatBeta_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="beta" select="count(FIELD[@name='MoffatBeta_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="betamax" select="count(FIELD[@name='MoffatBeta_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="residualsmin" select="count(FIELD[@name='Residuals_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="residuals" select="count(FIELD[@name='Residuals_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="residualsmax" select="count(FIELD[@name='Residuals_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pffwhmmin" select="count(FIELD[@name='FWHM_PixelFree_WCS_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pffwhm" select="count(FIELD[@name='FWHM_PixelFree_WCS_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pffwhmmax" select="count(FIELD[@name='FWHM_PixelFree_WCS_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfellipticitymin" select="count(FIELD[@name='Ellipticity_PixelFree_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfellipticity" select="count(FIELD[@name='Ellipticity_PixelFree_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfellipticitymax" select="count(FIELD[@name='Ellipticity_PixelFree_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfbetamin" select="count(FIELD[@name='MoffatBeta_PixelFree_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfbeta" select="count(FIELD[@name='MoffatBeta_PixelFree_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfbetamax" select="count(FIELD[@name='MoffatBeta_PixelFree_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfresidualsmin" select="count(FIELD[@name='Residuals_PixelFree_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfresiduals" select="count(FIELD[@name='Residuals_PixelFree_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pfresidualsmax" select="count(FIELD[@name='Residuals_PixelFree_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="asymmetrymin" select="count(FIELD[@name='Asymmetry_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="asymmetry" select="count(FIELD[@name='Asymmetry_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="asymmetrymax" select="count(FIELD[@name='Asymmetry_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="noiseqareamin" select="count(FIELD[@name='Area_Noise_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="noiseqarea" select="count(FIELD[@name='Area_Noise_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="noiseqareamax" select="count(FIELD[@name='Area_Noise_Max']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" onclick="showhideTable('psfext')" title="click to expand">
     PSF stats per Extension&nbsp;&darr;
    </BUTTON>
    <TABLE class="sortable" id="psfext" BORDER="2" style="display: none">
     <TR>
      <TH>Ext.</TH>
      <TH>N<sub>loaded</sub></TH>
      <TH>N<sub>accepted</sub></TH>
      <TH>Half-Light diam.</TH>
      <TH>Sampling</TH>
      <TH>&chi;<sup>2</sup>/d.o.f.</TH>
      <TH>FWHM</TH>
      <TH>FWHM<BR/>
       <elm>pixel free</elm>
      </TH>
      <TH>Ellipticity</TH>
      <TH>Ellipticity<BR/>
       <elm>pixel free</elm>
      </TH>
      <TH style="white-space: nowrap">Moffat &Beta;</TH>
      <TH style="white-space: nowrap">Moffat &Beta;<BR/>
       <elm>pixel free</elm>
      </TH>
      <TH>Residuals</TH>
      <TH>Residuals<BR/>
       <elm>pixel free</elm>
      </TH>
      <TH>Asymmetry</TH>
      <TH>Noise Equivalent Area</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
<!-- Extension number -->
        <td align="center">
         <el><xsl:value-of select="TD[$ext]"/></el>
        </td>
<!-- Number of sources loaded -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$nload],'######0.0')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="TD[$nloadmin]"/></elmin>
	  - <elmax><xsl:value-of select="TD[$nloadmax]"/></elmax>
         </elm>
        </td>
<!-- Number of sources accepted -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$naccep],'######0.0')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="TD[$naccepmin]"/></elmin>
	  - <elmax><xsl:value-of select="TD[$naccepmax]"/></elmax>
         </elm>
        </td>
<!-- Half-light diameter -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$fluxrad],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$fluxradmin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$fluxradmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Sampling -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$sampling],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$samplingmin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$samplingmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Chi2 --> 
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$chi2], '######0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$chi2min],'######0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$chi2max],'######0.00')"/></elmax>
         </elm>
        </td>
<!-- FWHM -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$fwhm],'##0.00')"/>"</el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$fwhmmin],'##0.00')"/>"</elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$fwhmmax],'##0.00')"/>"</elmax>
         </elm>
        </td>
<!-- FWHM pixel-free -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pffwhm],'##0.00')"/>"</el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pffwhmmin],'##0.00')"/>"</elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$pffwhmmax],'##0.00')"/>"</elmax>
         </elm>
        </td>
<!-- Ellipticity -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$ellipticity],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$ellipticitymin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$ellipticitymax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Ellipticity pixel-free -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pfellipticity],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pfellipticitymin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$pfellipticitymax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Beta -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$beta],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$betamin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$betamax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Beta pixel-free -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pfbeta],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pfbetamin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$pfbetamax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Residuals -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$residuals],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$residualsmin],'##0.00')"/></elmin>
 	  - <elmax><xsl:value-of select="format-number(TD[$residualsmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Residuals pixel-free -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$pfresiduals],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$pfresidualsmin],'##0.00')"/></elmin>
 	  - <elmax><xsl:value-of select="format-number(TD[$pfresidualsmax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Asymmetry -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$asymmetry],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$asymmetrymin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$asymmetrymax],'##0.00')"/></elmax>
         </elm>
        </td>
<!-- Noise Equivalent Area -->
        <td align="center">
         <el><xsl:value-of select="format-number(TD[$noiseqarea],'##0.00')"/></el>
         <br />
         <elm>
          <elmin><xsl:value-of select="format-number(TD[$noiseqareamin],'##0.00')"/></elmin>
	  - <elmax><xsl:value-of select="format-number(TD[$noiseqareamax],'##0.00')"/></elmax>
         </elm>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

<!-- ********************** XSL template for Config File ********************** -->
  <xsl:template name="Config">
   <p>
    <BUTTON type="button" onclick="showhideTable('config')" title="click to expand">
     Configuration File: <B><xsl:value-of select="PARAM[@name='Prefs_Name']/@value"/></B>&nbsp;&darr;
    </BUTTON>
    <TABLE id="config" class="sortable" style="display: none">
     <TR>
      <TH>Config Parameter</TH>
      <TH>Value</TH>
     </TR>
     <xsl:for-each select="PARAM[position()>2]">
      <tr>
       <td><el><xsl:value-of select="@name"/></el></td>
       <td><el><xsl:value-of select="@value"/></el></td>
      </tr>
     </xsl:for-each>
    </TABLE>
   </p>
   <p>
    <BUTTON type="button" onclick="showhideTable('commandline')" title="click to expand">
     Command Line&nbsp;&darr;
    </BUTTON>
    <TABLE id="commandline" style="display: none">
     <TR>
      <TD style="font-size: 80%;"><el><xsl:value-of select="PARAM[@name='Command_Line']/@value"/></el></TD>
     </TR>
    </TABLE>
   </p>
  </xsl:template>

<!-- ********************** XSL template for Warnings ********************** -->
  <xsl:template name="Warnings">
   <xsl:variable name="date" select="count(FIELD[@name='Date']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="time" select="count(FIELD[@name='Time']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="msg" select="count(FIELD[@name='Msg']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" onclick="showhideTable('warnings')" title="click to expand">
     Warnings (limited to the last 1000)&nbsp;&darr;
    </BUTTON>
    <TABLE id="warnings" class="sortable" style="display: none">
     <TR style="font-size: 80%;">
      <TH>Date</TH>
      <TH>Time</TH>
      <TH>Message</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td >
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td>
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$msg]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

 <xsl:template name="Rest">
</xsl:template>
</xsl:stylesheet>

