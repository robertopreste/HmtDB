
function populateCountries(s1, s2) {
    // Inserts countries values in the country list based on the selected continent
    s1 = document.getElementById(s1);
    s2 = document.getElementById(s2);
    var optionArray;

    s2.innerHTML = "--All Countries--";
    
    if (s1.value == "AF") {
        optionArray = ["|--All Countries--", 
        "6|Algeria", 
        "7|Angola", 
        "15|Barbados", 
        "22|Botswana", 
        "27|Burkina Faso", 
        "31|Cameroon", 
        "34|Central African Republic", 
        "37|Chad", 
        "43|Congo", 
        "58|Egypt", 
        "61|Eritrea", 
        "63|Ethiopia", 
        "69|Gabon", 
        "70|Gambia", 
        "73|Ghana", 
        "77|Guinea-Bissau", 
        "78|Haiti", 
        "97|Kenya", 
        "106|Libya", 
        "110|Madagascar", 
        "114|Mauritania", 
        "115|Mauritius", 
        "124|Morocco", 
        "125|Mozambique", 
        "127|Namibia", 
        "139|Niger", 
        "140|Nigeria", 
        "172|Sao Tome and Principe", 
        "176|Senegal", 
        "179|Sierra Leone", 
        "183|Somalia", 
        "184|South Africa", 
        "190|Sudan", 
        "196|United Republic of Tanzania", 
        "202|Tunisia", 
        "205|Uganda", 
        "207|Undefined African Country", 
        "228|Western Sahara", 
        "231|Zambia", 
        "232|Zimbabwe", 
        "244|Benin", 
        "248|British Indian Ocean Territory", 
        "249|Burundi", 
        "250|Cabo Verde", 
        "254|Comoros", 
        "256|Democratic Republic of the Congo", 
        "257|Cote d'Ivoire", 
        "260|Djibouti", 
        "262|Equatorial Guinea", 
        "265|French Southern Territories", 
        "272|Guinea", 
        "285|Lesotho", 
        "286|Liberia", 
        "291|Malawi", 
        "293|Mali", 
        "297|Mayotte", 
        "310|Reunion", 
        "311|Rwanda", 
        "313|Saint Helena Ascension and Tristan da Cunha", 
        "320|Seychelles", 
        "323|South Sudan", 
        "326|Swaziland", ];
    } else if (s1.value == "AM") {
        optionArray = ["|--All Countries--", 
        "4|Alabama", 
        "8|Argentina", 
        "13|Bahamas", 
        "18|Bermuda", 
        "20|Plurinational State of Bolivia", 
        "23|Brazil", 
        "29|California", 
        "32|Canada", 
        "33|Cayman Islands", 
        "39|Chile", 
        "41|Colombia", 
        "45|Costa Rica", 
        "48|Cuba", 
        "52|Delaware", 
        "54|Dominican Republic", 
        "57|Ecuador", 
        "59|El Salvador", 
        "67|Florida", 
        "75|Greenland", 
        "76|Guatemala", 
        "79|Hawai", 
        "83|Illinois", 
        "93|Kansas", 
        "96|Kentuky", 
        "113|Massachusetts", 
        "116|Mexico", 
        "117|Michigan", 
        "119|Minnesota", 
        "120|Missouri", 
        "123|Montana", 
        "128|Nebraska", 
        "132|New Hampshire", 
        "133|New Jersey", 
        "134|New Mexico", 
        "135|New York", 
        "138|Nicaragua", 
        "142|North Carolina", 
        "143|North Dakota", 
        "148|Ohio", 
        "149|Oklahoma", 
        "153|Panama", 
        "155|Paraguay", 
        "156|Pennsylvania", 
        "157|Peru", 
        "162|Puerto Rico", 
        "167|Rhode Island", 
        "185|South Carolina", 
        "186|South Dakota", 
        "197|Tennessee", 
        "198|Texas", 
        "208|Undefined American Country", 
        "215|Uruguay", 
        "216|United States of America", 
        "219|Bolivian Republic of Venezuela", 
        "220|Vermont", 
        "223|Virginia", 
        "226|Washington", 
        "229|Wisconsin", 
        "234|Alaska", 
        "237|Anguilla", 
        "238|Antigua and Barbuda", 
        "239|Arizona", 
        "240|Arkansas", 
        "241|Aruba", 
        "243|Belize", 
        "246|Bonaire Sint Eustatius and Saba", 
        "247|Bouvet Island", 
        "253|Colorado", 
        "255|Connecticut", 
        "258|Curacao", 
        "259|District of Columbia", 
        "261|Dominica", 
        "263|Falkland Islands", 
        "264|French Guiana", 
        "266|Georgia (US)", 
        "268|Grenada", 
        "269|Guadeloupe", 
        "273|Guyana", 
        "276|Honduras", 
        "278|Idaho", 
        "279|Indiana", 
        "280|Iowa", 
        "282|Jamaica", 
        "288|Louisiana", 
        "290|Maine", 
        "295|Martinique", 
        "296|Maryland", 
        "298|Mississippi", 
        "301|Montserrat", 
        "303|Nevada", 
        "307|Oregon", 
        "312|Saint Barthelemy", 
        "314|Saint Kitts and Nevis", 
        "315|Saint Lucia", 
        "316|Saint Martin", 
        "317|Saint Pierre and Miquelon", 
        "318|Saint Vincent and the Grenadines", 
        "321|Sint Marteen", 
        "322|South Georgia and the South Sandwich Islands", 
        "324|Suriname", 
        "330|Trinidad and Tobago", 
        "332|Turks and Caicos Islands", 
        "334|Utah", 
        "335|Virgin Islands (British)", 
        "336|Virgin Islands (US)", 
        "337|West Virginia", 
        "338|Wyoming", ];
    } else if (s1.value == "AS") {
        optionArray = ["|--All Countries--", 
        "2|Afghanistan", 
        "9|Armenia", 
        "12|Azerbaijan", 
        "14|Bangladesh", 
        "25|Brunei Darussalam", 
        "30|Cambodia", 
        "40|China", 
        "56|Timor-Leste", 
        "71|Georgia", 
        "84|India", 
        "85|Indonesia", 
        "86|Islamic Republic of Iran", 
        "87|Iraq", 
        "89|Israel", 
        "91|Japan", 
        "92|Jordan", 
        "95|Kazakhstan", 
        "99|Korea", 
        "101|Kuwait", 
        "102|Kyrgyzstan", 
        "103|Laos People's Democratic Republic", 
        "105|Lebanon", 
        "111|Malaysia", 
        "122|Mongolia", 
        "126|Myanmar", 
        "129|Nepal", 
        "150|Oman", 
        "151|Pakistan", 
        "152|State of Palestine", 
        "158|Philippines", 
        "163|Qatar", 
        "165|Singapore", 
        "173|Saudi Arabia", 
        "188|Sri Lanka", 
        "194|Syrian Arab Republic", 
        "195|Taiwan", 
        "199|Thailand", 
        "203|Turkey", 
        "209|Undefined Asian Country", 
        "213|United Arab Emirates", 
        "217|Uzbekistan", 
        "221|Viet Nam", 
        "230|Yemen", 
        "242|Bahrain", 
        "245|Bhutan", 
        "277|Hong Kong", 
        "284|Democratic People's Republic of Korea", 
        "289|Macao", 
        "292|Maldives", 
        "327|Tajikistan", 
        "328|Togo", 
        "331|Turkmenistan", ];
    } else if (s1.value == "EU") {
        optionArray = ["|--All Countries--", 
        "5|Albania", 
        "11|Austria", 
        "16|Belarus", 
        "17|Belgium", 
        "21|Bosnia and Herzegovina", 
        "26|Bulgaria", 
        "47|Croatia", 
        "49|Cyprus", 
        "50|Czechia", 
        "53|Denmark", 
        "62|Estonia", 
        "64|Faroe Islands", 
        "66|Finland", 
        "68|France", 
        "72|Germany", 
        "74|Greece", 
        "80|Hungary", 
        "82|Iceland", 
        "88|Ireland", 
        "90|Italy", 
        "104|Latvia", 
        "107|Lithuania", 
        "108|Luxembourg", 
        "109|Macedonia", 
        "112|Malta", 
        "121|Republic of Moldova", 
        "130|Netherlands", 
        "147|Norway", 
        "159|Poland", 
        "161|Portugal", 
        "168|Romania", 
        "170|Russian Federation", 
        "177|Serbia", 
        "180|Slovakia", 
        "181|Slovenia", 
        "187|Spain", 
        "191|Sweden", 
        "192|Switzerland", 
        "206|Ukraine", 
        "211|Undefined European Country", 
        "214|United Kingdom of Great Britain and Northern Ireland", 
        "233|Aland Islands", 
        "236|Andorra", 
        "267|Gibraltar", 
        "271|Guernsey", 
        "275|Holy See", 
        "281|Isle of Man", 
        "283|Jersey", 
        "287|Liechtenstein", 
        "299|Monaco", 
        "300|Montenegro", 
        "319|San Marino", 
        "325|Svalbard and Jan Mayen", ];
    } else if (s1.value == "OC") {
        optionArray = ["|--All Countries--", 
        "10|Australia", 
        "44|Cook Islands", 
        "65|Fiji", 
        "98|Kiribati", 
        "118|Federal States of Micronesia", 
        "136|New Zealand", 
        "141|Niue", 
        "154|Papua New Guinea", 
        "160|French Polynesia", 
        "171|Samoa", 
        "182|Solomon Islands", 
        "201|Tonga", 
        "204|Tuvalu", 
        "212|Undefined Oceanian Country", 
        "218|Vanuatu", 
        "225|Wallis and Futuna", 
        "235|American Samoa", 
        "251|Christmas Island", 
        "252|Cocos (Keeling) Islands", 
        "270|Guam", 
        "274|Heard Island and McDonald Islands", 
        "294|Marshall Islands", 
        "302|Nauru", 
        "304|New Caledonia", 
        "305|Norfold Island", 
        "306|Northern Mariana Islands", 
        "308|Palau", 
        "309|Pitcairn", 
        "329|Tokelau", 
        "333|United States Minor Outlying Islands", ]; 
    }

    for (var option in optionArray) {
        var pair = optionArray[option].split("|");
        var newOption = document.createElement("option");

        newOption.value = pair[0];
        newOption.innerHTML = pair[1];
        s2.options.add(newOption);
    }
}



function openCenteredPopup(target, w, h) {
    // Opens centered popups for elements like Statistics, People, ecc.
    var t = Math.floor((screen.height-h)/2);
    var l = Math.floor((screen.width-w)/2);
    var stile = "top=" + t + ", left=" + l + ", width=" + w + ", height=" + h + ", status=no, menubar=no, toolbar=no, scrollbar=no";
    window.open(target, "", stile);
}

