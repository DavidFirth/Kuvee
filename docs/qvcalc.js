/*
quasiVariance(covmat): Calculates quasi variances from a covariance matrix
Throws invalidConVar
*/
quasiVariance = function (covmat) {
    var covStyle = matrixType(covmat);

    if(covStyle == ""){
        throw new Error("The matrix is not of one of the 3 viable formats.");
	   validFlag = 0;
    }
    if(covStyle == "Errors and correlation matrix"){
        covmat = corrToCov(covmat);
    }
    if(!checkCorr(covmat)){
            throw new Error("The matrix implies correlations not in [-1,1],cannot be a valid covariance/correlation matrix.");
    }

    var xmat = desMat(covmat.length);
    var vvec = simConVar(covmat);
    //Debug information
    //console.log("SC Variances:" + vvec);
    //console.log("SC log Variances:" + vvec.map(Math.log));
    //document.getElementById("debugOut").innerHTML="<h2>(Debug information:)</h2><br>contrast variance:" + vecApply(vvec, Math.log) + "<br>simCon:" + simCon(covmat.length) + "<br><br>desMat:<br>" + matrixToHTML(xmat);

    //Observed values are the log contrast variances
    if(vvec.isPositive()) var yvec = vvec.map(Math.log);
    else {
        throw new Error("Some simple contrast variances were non-positive. Are you sure you have provided a valid covariance/correlation matrix?");
    }
    //Initialise glm
    var glm_model = GLM(GLM.families.Gaussian(GLM.links.linkBuilder(Math.exp,Math.log,Math.exp)));
    //Fit glm
    glm_model.fit(yvec, xmat);
    var resid = yvec.minus(glm_model.predict(xmat));
    //console.log("Resids: " + resid);
    var relerrs = resid.map(Math.exp);
    relerrs = relerrs.map(Math.sqrt);
    relerrs = relerrs.map(function(x){return 1-x;});
    //console.log("Relerrs:" + relerrs.toFixed(5));

    return {SE: getDiag(covmat).map(Math.sqrt),
            quasiVar: glm_model.weights,
            quasiSE: glm_model.weights.map(Math.sqrt),
            relerrs: relerrs,
            covmat: symmetrify(covmat)};
}

worstErrors = function (array, qvar){
    var covmat = new Matrix(array);
    var levels = covmat.getRowDimension();
    var J = new Matrix(levels, levels, 1); //levels x levels matrix of 1s
    var ones = rep(1, levels);
    var t1 = outerProduct(array[0], ones);
    var t2 = outerProduct(ones, array[0]);
    var r_covmat = covmat.plus(J.timesScalar(array[0][0])).minus(t1).minus(t2)
    var r_covmat = r_covmat.getMatrix(1, levels-1, 1, levels-1);
    var qvmat = diag(qvar);
    var r_qvmat = qvmat.plus(J.timesScalar(qvar[0])).getMatrix(1, levels-1, 1, levels-1);
    var L = r_covmat.chol().getL();
    var L_inv = L.inverse();
    var emat = L_inv.times(r_qvmat).times(L_inv.transpose());
    var evals = emat.eig().getD();
    var minqv = Math.sqrt(evals.getDiag().min())-1;
    var maxqv = Math.sqrt(evals.getDiag().max())-1;
    return {min: minqv, max:maxqv};
}

outerProduct = function (vec1,vec2) {
    var n = vec1.length;
    var m = vec2.length;
    var mat = new Matrix(n,m);
    for (var i = 0; i < n; i++) {
        for (var j = 0; j < m; j++) {
            mat.set(i,j, vec1[i]*vec2[j])
        }
    }
    return mat;
}

//rep: a vector length n with each entry being "val"
rep = function (val, n) {
   var arr = new Array(n);
   for(i=0; i<n; i++) arr[i] = val;
   return arr;
}

//diag: Make a diagonal matrix from a vector
diag = function (vec) {
   n = vec.length;
   var mat = new Matrix(n,n,0);
   for(i=0; i<n; i++) mat.set(i,i, vec[i]);
   return mat;
}



//symmetrify: makes a triangle of cov matrix into square
symmetrify = function (arr) {
    var m = arr.length;
    var symm = arr.slice();
    for(var i = 1; i < m; i++){
    for(var j = 0; j < i; j++){
        symm[j].push(arr[i][j]);
    }}
   return symm;
}

//getDiag: Get the diagonal of a 2d matrix
getDiag = function (arr) {
    var n = arr.length;
    var vec = new Array(n);
    for(i=0; i<n; i++) vec[i] = arr[i][i];
    return vec;
}



//matrixType: Determines whether the matrix is of 1 of the 3 allowed formats
matrixType = function (mat) {
    var type = "";
    var m = mat.length;
    var n = new Array(m);
    var squareFlag = 1;
    var triangleFlag = 1;
    var vecTriFlag = 1;
    if(mat[0].length != m) vecTriFlag = 0;
    for (var i = 0; i < m; i++){
        n[i] = mat[i].length;
        if(n[i] != m) squareFlag = 0;
        if(n[i] != i+1) triangleFlag = 0;
        if(n[i] != i && i != 0) vecTriFlag = 0;
    }
    if(squareFlag == 1) type="Square covariance matrix";
    if(triangleFlag == 1) type="Lower covariance matrix";
    if(vecTriFlag == 1) type="Errors and correlation matrix";
    return type;
}

//corrToCov: if the data given is of the "Errors and correlation matrix" form,
//            this function converts it to the "lower covariance matrix" form
corrToCov = function (mat) {
    var m = mat.length;
    var Rmat = new Array(m);
    for(var i = 1; i < m; i++){
        Rmat[i] = mat[i];
        Rmat[i].push(1);
    }
    Rmat[0] = new Array(1);
    Rmat[0][0] = mat[0][0];
    for(var i = 0; i < m; i++){
    for(var j = 0; j <= i; j++){
        Rmat[i][j] = Rmat[i][j]*mat[0][i]*mat[0][j];
    }}
    return Rmat;
}


//checkCorr: checks whether a given covariance matrix has invalid correlation coefficients.
checkCorr = function (mat) {
    var flag = 1;
    var m = mat.length;
    for(var i = 1; i < m; i++){
    for(var j = 0; j < i; j++){
        if(mat[i][i] == 0 || mat[j][j] == 0 || mat[i][j] == 0) var corr = 0;
        else var corr = mat[i][j]/(Math.sqrt(mat[i][i])*Math.sqrt(mat[j][j]));
        //console.log("Corr[" + i + "][" + j + "]: " + corr);
        if (corr > 1 || corr < -1){
            flag = 0;
        }
    }}
    return flag;
}

//simConVar: Finds the variance of all simple contrasts given the covariance matrix (only needs lower triangle)
simConVar = function (cov) {
    var n = cov.length;
    var con = simCon(n);
    var variance = new Array(binom(n,2));
    for (var i=0; i < binom(n,2); i++){
        j = con[i][0];
        k = con[i][1];
        var tmp = Number(cov[j][j])+Number(cov[k][k])-Number(2*cov[k][j]);
        variance[i] = tmp;
    }
    return variance;
}

//binom: calculates "n choose k"
binom = function (n,k) {
    var coeff = 1;
    for (var i = n-k+1; i <= n; i++) coeff *= i;
    for (var i = 1;     i <= k; i++) coeff /= i;
    return coeff;
}

//simCon: generates simple contrasts in numbers
simCon = function (levels) {
    var k = 0;
    var mat = new Array(binom(levels,2));
    for (var i = 0; i < levels - 1; i++){
        for (var j = i + 1; j < levels; j++){
            mat[k] = [i,j];
            k++;
        }
    }
    return mat;
}

//desMat: generates simple contrasts in binary
desMat = function (levels) {
    var k = 0;
    var tmp = new Array(levels);
    var des = new Array(binom(levels,2));
    for (var i = 0; i < levels; i++){
        tmp[i] = 0;
    }
    for (var i = 0; i < levels - 1; i++){
        for (var j = i + 1; j < levels; j++){
            des[k] = tmp.slice();
            des[k][i] = 1;
            des[k][j] = 1;
            k++;
        }
    }
    return des;
}

/******************************/
/****** Array functions *******/
/******************************/

Array.prototype.max = function () {
  return Math.max.apply(null, this)
}

Array.prototype.min = function () {
  return Math.min.apply(null, this)
}

//minus: Takes "arr" away from this array element-wise
Array.prototype.minus = function (arr) {
    if (this.length != arr.length) throw "nEqErr";
    var vec = new Array(this.length);
    for(var i = 0; i < this.length; i++){
       vec[i] = Number(this[i]) - Number(arr[i]);
    }
    return vec;
}

//plus: adds "arr" from this array element-wise
Array.prototype.plus = function (arr) {
    if (this.length != arr.length) throw "nEqErr";
    var vec = new Array(this.length);
    for(var i = 0; i < this.length; i++){
       vec[i] = Number(this[i]) + Number(arr[i]);
    }
    return vec;
}

//plus: adds "arr" from this array element-wise
Array.prototype.times = function (arr) {
    if (this.length != arr.length) throw "nEqErr";
    var vec = new Array(this.length);
    for(var i = 0; i < this.length; i++){
       vec[i] = Number(this[i])*Number(arr[i]);
    }
    return vec;
}

//isPositive: checks if an array has only positive elements
Array.prototype.isPositive = function () {
    var tmp = 1;
    for(var i=0; i < this.length; i++){
        if(this[i] <= 0) tmp = 0;
    }
    return tmp;
}

Array.prototype.toFixed = function (n) {
    var k = this.length;
    var i = 0;
    var Rvec = new Array(k);
    while(i < k) {
            Rvec[i] = this[i].toFixed(n);
            i++;
    }
    return Rvec;
}

Array.prototype.euclidNorm = function () {
    var sum = 0;
    for(var i = 0; i < this.length; i++){
        sum += (this[i])*(this[i])
    }
    return Math.sqrt(sum);
}

Array.prototype.penalty = function (v) {
    var penalty = new Array(v.length);
    var k = this.length;
    var i = 0, j = 1, l=0;
    while ( i < k-1 ){ //"While" loops faster than "for" loops in JS
        j = i+1;
        while ( j < k ){
            penalty[l] = Math.pow(Math.log((this[i] + this[j])/v[l]),2);
            j++;
            l++;
        }
        i++;
    }
    return penalty;
}

Array.prototype.simCon = function () {
    var k = this.length;
    var sc = new Array(binom(k,2));
    var i = 0, j = 1, l=0;
    while ( i < k-1 ){ //"While" loops faster than "for" loops in JS
        j = i+1;
        while ( j < k ){
            sc[l] = this[i] + this[j];
            j++;
            l++;
        }
        i++;
    }
    return sc;
}

Array.prototype.scaleBy = function (alpha) {
    var result = new Array(this.length);
    for(var i=0; i < this.length; i++){
        result[i] = this[i]*alpha;
    }
    return result;
}

Array.prototype.residuals = function (v) {
    var penalty = new Array(v.length);
    var k = this.length;
    var i = 0, j = 1, l=0;
    while ( i < k-1 ){ //"While" loops faster than "for" loops in JS
        j = i+1;
        while ( j < k ){
            penalty[l] = (this[i] + this[j])/v[l] - 1;
            j++;
            l++;
        }
        i++;
    }
    return penalty;
}

Array.prototype.penalties = function (v) {
    var penalty = new Array(v.length);
    var k = this.length;
    var i = 0, j = 1, l=0;
    while ( i < k-1 ){ //"While" loops faster than "for" loops in JS
        j = i+1;
        while ( j < k ){
            penalty[l] = Math.pow(Math.log((this[i] + this[j])/v[l]),2);
            j++;
            l++;
        }
        i++;
    }
    return penalty;
}

Array.prototype.divideBy = function (v) {
    var k = this.length, i = 0;
    if(k != v.length) throw new Error ("DivideBy: Vectors of different length - cannot do elementwise division.");
    var result = new Array(k);
    while (i < k){
        if(v[i] == 0) throw new Error ("DivideBy: Divisor is 0 - division is undefined.")
        result[i] = this[i] / v[i];
        i++;
    }
    return result;
}
/*******************************/
/****** Matrix functions *******/
/*******************************/

Array.prototype.isNull = function (tol) {
    if(this.map(Math.abs).max() < 1/Math.pow(10,tol)) return 1;
    else return 0;
}

Matrix.prototype.getDiag = function () {
    var n = this.getRowDimension();
    var vec = new Array(n);
    for(i=0; i<n; i++) vec[i] = this.get(i,i);
    return vec;
}

//Matrix.toTable(): Method to write to write a Matrix to a HTML table format
Matrix.prototype.toTable = function () {
    var str = "<table>";
    for (var i = 0; i < this.getRowDimension(); i++){
        str = str + "<tr>";
        for (var j = 0; j < this.getColumnDimension(); j++){
            str = str + "<td>  " + this.get(i,j) + "  </td>";
        }
        str = str + "</tr>";
    }
    return str+ "</table></br>";
}



/* An Example:
0 0 0 0 0
0 0.0533 0.0428 0.0390 0.0404
0 0.0428 0.1824 0.0384 0.0412
0 0.0390 0.0384 0.1428 0.0387
0 0.0404 0.0412 0.0387 0.0941 

should give

                SE quasiSE quasiVar
      level 1 0.000   0.201   0.0404
      level 2 0.231   0.113   0.0127
      level 3 0.427   0.374   0.1401
      level 4 0.378   0.324   0.1049
      level 5 0.307   0.232   0.0539


This currently crashes the glm:
0 0 0 0 0 0 0 0 0 0
0 1 0 0 0 0 0 0 3 0
0 0 2 0 0 0 0 1 0 0
0 0 0 3 0 0 0 0 0 0
0 0 0 0 4 0 0 0 0 1
0 0 0 0 0 5 0 0 0 0
0 0 0 0 0 0 6 0 0 0
0 0 1 0 0 0 0 7 0 0
0 3 0 0 0 0 0 0 8 0
0 0 0 0 1 0 0 0 0 9

Replacing the 3 with something <= 2.413959 works

0
0 1
0 0 2
0 0 0 3
0 0 0 0 4
0 0 0 0 0 5
0 0 0 0 0 0 6
0 0 1 0 0 0 0 7
0 2.413959 0 0 0 0 0 0 8
0 0 0 0 1 0 0 0 0 9

the same effect can be seen with

0
0 1
0 1.99 3

replacing the 1.99 with a 1.5 solves the issue.

Happens in linearsolve on second loop of IRLS
*/
