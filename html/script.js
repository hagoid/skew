let checkboxAuto
let checkboxCoset
let checkboxOther
let autoVisible
let cosetVisible
let otherVisible
let modulo
let moduloInput
let selected

$.fn.countSelected = function() {
    return this.each(function() {
        let elem = $(this)
        elem.find('.count').text(elem.find('.skew:visible').length)
    })
}

$.fn.represent = function() {
    return this.each(function() {
        let elem = $(this)
        let data = elem.data('repr')
        let repr = '(' + data.map(e => e.map(f => f % modulo).join(', ')).join(')(') + ')'
        if (repr === '()' || repr === '(0)' || repr === '(1)') repr = 'Id(sn)'
        elem.find('.repr').text('with representant ' + repr)
    })
}

$.fn.linking = function() {
    let params = {
        'auto': autoVisible,
        'coset': cosetVisible,
        'other': otherVisible,
        'mod': modulo,
    }
    return this.each(function() {
        let elem = $(this)
        let id = elem.prop('id')
        elem.find('a.link').prop('href', '?' + $.param(params) + '#' + id)
    })
}

let update = () => {
    $('.auto').toggle(autoVisible)
    $('.coset').toggle(cosetVisible)
    $('.other').toggle(otherVisible)
    $('.n').countSelected()
    $('.skew').represent().linking()
    if (selected) selected.toggleClass('selected', false)
    selected = $(window.location.hash + ",." + window.location.hash.split('#')[1])
    selected.toggleClass('selected', true)
}

$(document).ready(() => {
    checkboxAuto = $('#checkbox-auto')
    checkboxCoset = $('#checkbox-coset')
    checkboxOther = $('#checkbox-other')
    moduloInput = $('#modulo')
    let urlParams = new URLSearchParams(window.location.search)
    if (urlParams.has('auto')) checkboxAuto.prop('checked', urlParams.get('auto') === 'true')
    if (urlParams.has('coset')) checkboxCoset.prop('checked', urlParams.get('coset') === 'true')
    if (urlParams.has('other')) checkboxOther.prop('checked', urlParams.get('other') === 'true')
    if (urlParams.has('mod')) moduloInput.prop('value', urlParams.get('mod'))
    autoVisible = checkboxAuto.prop('checked')
    cosetVisible = checkboxCoset.prop('checked')
    otherVisible = checkboxOther.prop('checked')
    modulo = moduloInput.prop('value')
    checkboxAuto.change(() => {
        autoVisible = checkboxAuto.prop('checked')
        update()
    })
    checkboxCoset.change(() => {
        cosetVisible = checkboxCoset.prop('checked')
        update()
    })
    checkboxOther.change(() => {
        otherVisible = checkboxOther.prop('checked')
        update()
    })
    moduloInput.change(() => {
        modulo = moduloInput.prop('value')
        update()
    })
    $( window ).on( 'hashchange', function( e ) {
        if (selected) selected.toggleClass('selected', false)
        selected = $(window.location.hash)
        selected.toggleClass('selected', true)
    } );
    update()
    $(document).keydown(function(e){
        if (e.which == 37) {
            $('a#prev')[0].click()
        } else if (e.which == 39) {
            $('a#next')[0].click()
        }
    })
    window.onhashchange = function() {
        update()
    }
})

const cyrb53 = (vector, seed = 0) => {
    let h1 = 0xdeadbeef ^ seed,
        h2 = 0x41c6ce57 ^ seed;
    for (element of vector) {
        h1 = Math.imul(h1 ^ element, 2654435761);
        h2 = Math.imul(h2 ^ element, 1597334677);
    }

    h1 = Math.imul(h1 ^ (h1 >>> 16), 2246822507) ^ Math.imul(h2 ^ (h2 >>> 13), 3266489909);
    h2 = Math.imul(h2 ^ (h2 >>> 16), 2246822507) ^ Math.imul(h1 ^ (h1 >>> 13), 3266489909);

    return 4294967296 * (2097151 & h2) + (h1 >>> 0);
}