 /**
   * @return {string} The tooltip text
   */
  this.getTooltipText = function(last_parent, caller) {
    if (!last_parent && caller.hid === "parent")
      return "" + nb_parents + " parents";
    else if(caller.hid === "parent" || caller.hid === "sibling")
      return "" + nb_siblings + " siblings";
    else if(caller.hid === "current" && nb_children==0)
      return "The current object";
    else if(caller.hid === "children" || caller.hid === "current" )
      return "" + nb_children + " children";

/*
    if(this.hid === "parent")
      return "" + nb_parents + " parents";
    else if(this.hid === "sibling")
      return "" + nb_siblings + " siblings";
    else if(this.hid === "children")
      return "" + nb_children + " children";
    else if(this.hid === "current")
      return "The current object";
*/
  };

var SIMBAD = {};

/**
 * Definition d'une classe représentant un cercle
 * @param {int} xc     x du centre
 * @param {int} yc     y du centre
 * @param {int} radius Rayon du cercle
 * @param {int} hid  parent; sibling; children; current
 */
SIMBAD.Circle = function(xc, yc, radius, hid) {
  this.xc = xc;
  this.yc = yc;
  this.radius = radius;
  this.hid = hid;

  /**
   * Est ce que le point passé en parametre est contenu dans le cercle ?
   * @param  {int} x     Abscisse du point à tester
   * @param  {int} y     Ordonnée du point à tester
   * @param  int} scale  Permet de mettre à l'echelle le rayon en retranchant ou ajoutant une somme
   * @return {boolean}       Vrai ou faux
   */
  this.contains = function(x, y, scale) {
    var rad = this.radius + scale;
    square_dist = Math.pow((this.xc - x), 2) + Math.pow((this.yc - y), 2);
    return (square_dist <= Math.pow(rad, 2));
  };

  /**
   * @return {String} L'url à suivre lors d'un clic sur le cercle
   */
  this.getTargetUrl = function(target) {
    if(target === "parent") {
	return "$('input[value=\"parents\"]').click();";
	// URL + NbIdent=query_hlinks&submit=parents
    //  return "http://www.google.fr";
    }
    else if(target === "sibling") {
	return "$('input[value=\"siblings\"]').click();"
     // return "http://www.google.fr";
    }
    else if(target === "children" || target === "current"){
	return "$('input[value=\"children\"]').click();"
    }
      //return "http://www.google.fr";
  };


  /**
   * Fonction qui permet de dessiner le cercle
   */
  this.draw = function(last_parent) {
    last_parent = last_parent || false;

    var title = document.createElementNS("http://www.w3.org/2000/svg", "title");
    title.textContent = getTooltipText(last_parent, this);
    var shape = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    shape.setAttributeNS(null, "cx", this.xc);
    shape.setAttributeNS(null, "cy", this.yc);
    shape.setAttributeNS(null, "r", this.radius);

    var classe = "" + this.hid;
    if(last_parent)
      classe += " last_parent";
    shape.setAttributeNS(null, "class", classe); // classe CSS

//    if(!last_parent && this.hid !== "parent")
	//if (this.hid !== "parent")
	if (this.hid==="parent" && last_parent && nb_siblings>0)
      shape.setAttributeNS(null, "onclick", this.getTargetUrl("sibling") );
	else if (nb_children>0 && (this.hid==="current" ||this.hid==="children"))
      shape.setAttributeNS(null, "onclick", this.getTargetUrl("children")) ;
	else if (nb_siblings>0 && this.hid==="sibling")
      shape.setAttributeNS(null, "onclick", this.getTargetUrl("sibling") );
	else if (this.hid==="parent" && !last_parent && nb_siblings>0)
	shape.setAttributeNS(null, "onclick", this.getTargetUrl("parent") );


      //shape.setAttributeNS(null, "onclick", "location.href='" + this.getTargetUrl() + "'");

    shape.appendChild(title); // tooltip
    document.getElementById('hierarchy').appendChild(shape);
  };
}

