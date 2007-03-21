<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
	]>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- *********************** Global XSL template ************************** -->
 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <html>

<!-- HTML head  -->

   <head>

<!-- javascript -->

<!--  <script type="text/javascript" language="javascript"> -->

    <script src="http://terapix.iap.fr/cplt/xsl/sorttable.js"/>

    <style type="text/css">
     p.sansserif {font-family: sans-serif}
     body {background-color: white}
     mono {font-family: monospace}
     elen {font-family: monospace; font-size: 100%; font-weight: bold; color: green }
     elep {font-family: monospace; font-size: 100%; font-weight: bold; color: red }
     el {font-family: monospace; font-size: 100%; color: black}
     a {text-decoration: none}
     table.sortable a.sortheader
      {
      background-color:#FFEECC;
      color: black;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }
     table.sortable span.sortarrow
      {
      color: black;
      font-weight: bold;
      text-decoration: none;
      }
     table.sortable a.sortheader.sub
      {
      vertical-align: sub;
      }
     </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </head>

<!-- HTML body -->

    <BODY>
     <TABLE BORDER="0" CELLPADDING="0" CELLSPACING="0" WIDTH="100%">
      <TR>
       <TD ALIGN="LEFT">
        <TABLE BORDER="0">
         <TR>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixLogo.png" ALT="Terapix"/>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixTitle.png" ALT="Logo"/>
          </TD>
          <TD ALIGN="CENTER">
           <FONT color="#669933">
            <B> Processing summary</B>
           </FONT>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixPicture.gif" ALT="Terapix banner"/>
          </TD>
         </TR>
        </TABLE>
       </TD>
      </TR>
      <TR>
       <TD>
        <TABLE BORDER="0" WIDTH="100%" BGCOLOR="#000000">
         <TR>
          <TH BGCOLOR="#000000" ALIGN="LEFT"><FONT SIZE="-1" COLOR="#FFFFFF"> Home > Tools > Data reduction</FONT></TH>
         </TR>
        </TABLE>
       </TD>
      </TR>
     </TABLE>
    <xsl:call-template name="VOTable"/>
   </BODY>
  </html>
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
    <xsl:when test="@ID='SCAMP'">
     <xsl:call-template name="SCAMP"/>
    </xsl:when>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>

<!-- ********************** XSL template for SCAMP *********************** -->
 <xsl:template name="SCAMP">
  <xsl:for-each select="RESOURCE[@ID='MetaData']">
   <xsl:call-template name="RunInfo"/>
   <xsl:for-each select="TABLE[@ID='Fields']">
    <xsl:call-template name="Fields"/>
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
 </xsl:template>

<!-- ********************** XSL template for Fields ********************** -->
  <xsl:template name="Fields">
   <xsl:variable name="name" select="count(FIELD[@name='Catalog_Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ident" select="count(FIELD[@name='Image_Ident']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="next" select="count(FIELD[@name='NExtensions']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ndet" select="count(FIELD[@name='NDetect']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="headflag" select="count(FIELD[@name='Ext_Header']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="photflag" select="count(FIELD[@name='Photom_Flag']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="group" select="count(FIELD[@name='Group']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astrinstru" select="count(FIELD[@name='Astr_Instrum']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="photinstru" select="count(FIELD[@name='Phot_Instrum']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="coord" select="count(FIELD[@name='Field_Coordinates']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="radius" select="count(FIELD[@name='Max_Radius']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pixscale" select="count(FIELD[@name='Pixel_Scale']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ascont" select="count(FIELD[@name='AS_Contrast']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="xycont" select="count(FIELD[@name='XY_Contrast']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="zpcorr" select="count(FIELD[@name='ZeroPoint_Corr']/preceding-sibling::FIELD)+1"/>
   <p>
    <TABLE class="sortable" id="scamp" BORDER="2">
<!--
 <TR>
  <TH COLSPAN="3" BGCOLOR="#DDCC99" ALIGN="LEFT">
  <FONT SIZE="-1">
   <xsl:value-of select="ASTROM/nfield"/>
    input<xsl:if test="ASTROM/nfield &gt; 1">s</xsl:if>
  </FONT>
  </TH>
 </TR>
-->
     <TR>
      <TH BGCOLOR="#FFEECC">Filename</TH>
      <TH BGCOLOR="#FFEECC">Identifier</TH>
      <TH BGCOLOR="#FFEECC">Next</TH>
      <TH BGCOLOR="#FFEECC">Ndet</TH>
      <TH BGCOLOR="#FFEECC">Flags</TH>
      <TH BGCOLOR="#FFEECC">G</TH>
      <TH BGCOLOR="#FFEECC">A</TH>
      <TH BGCOLOR="#FFEECC">P</TH>
      <TH BGCOLOR="#FFEECC">alpha</TH>
      <TH BGCOLOR="#FFEECC">delta</TH>
      <TH BGCOLOR="#FFEECC">Radius</TH>
      <TH BGCOLOR="#FFEECC">Pixel scale</TH>
      <TH BGCOLOR="#FFEECC">A/S contrast</TH>
      <TH BGCOLOR="#FFEECC">X/Y contrast</TH>
      <TH BGCOLOR="#FFEECC">MagZP.corr</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ident]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$next]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ndet]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$headflag] = 'T'">
           <elen>H</elen>
          </xsl:when>
          <xsl:otherwise>
            <el>-</el>
          </xsl:otherwise>
         </xsl:choose>
         <xsl:choose>
          <xsl:when test="TD[$photflag] = 'T'">
            <elen>P</elen>
          </xsl:when>
          <xsl:otherwise>
            <el>-</el>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$group]"/></el>
        </td>
        <td align="left" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astrinstru]"/></el>
        </td>
        <td align="left" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$photinstru]"/></el>
        </td>
<!-- Alpha -->
        <td align="center" BGCOLOR="#EEEEEE">
         <el>
          <xsl:variable name="alpha" select="substring-before(string(TD[$coord]), ' ')"/>
          <xsl:value-of
		select='concat(format-number(floor($alpha div 15.0), "00"),":",
		format-number(floor(($alpha * 4) mod 60.0), "00"),":",
		format-number(floor(($alpha * 240.0) mod 60.0), "00.00"))'/>
         </el>
        </td>
<!-- Delta -->
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:variable name="delta" select="substring-after(string(TD[$coord]), ' ')"/>
         <el>
          <xsl:choose>
           <xsl:when test="$delta &lt; 0.0">
            <xsl:value-of
		select='concat("-", format-number(floor(-$delta), "00"),":",
		format-number(floor((-$delta * 60) mod 60.0), "00"),":",
		format-number(floor((-$delta * 3600.0) mod 60.0), "00.0"))'/>
           </xsl:when>
           <xsl:otherwise>
            <xsl:value-of
		select='concat("+", format-number(floor($delta), "00"),":",
		format-number(floor(($delta * 60) mod 60.0), "00"),":",
		format-number(floor(($delta * 3600.0) mod 60.0), "00.0"))'/>
           </xsl:otherwise>
          </xsl:choose>
         </el>
        </td>
<!-- Radius -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$radius],'##0.000')"/>&amin;</el>
        </td>
<!-- Pixel scale --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:variable name="pix1" select="number(substring-before(string(TD[$pixscale]), ' '))"/>
         <xsl:variable name="pix2" select="number(substring-after(string(TD[$pixscale]), ' '))"/>
         <el><xsl:value-of select="format-number(($pix1+$pix2) div 2.0, '##0.0000')"/>&asec;</el>
        </td>
<!-- A/S contrast --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$ascont], '##0.0')"/></el>
        </td>
<!-- X/Y contrast --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$xycont] &lt; 2.0">
           <elep><xsl:value-of select="format-number(TD[$xycont], '##0.0')"/></elep>
          </xsl:when>
          <xsl:otherwise>
            <elen><xsl:value-of select="format-number(TD[$xycont], '##0.0')"/></elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
<!-- Zero-Point correction --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="contains(TD[$zpcorr],'e')">
           <el><xsl:value-of select="format-number('0', '#0.000')"/></el>
          </xsl:when>
          <xsl:otherwise>
           <el><xsl:value-of select="format-number(TD[$zpcorr], '#0.000')"/></el>
          </xsl:otherwise>
         </xsl:choose>
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
