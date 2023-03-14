let checkboxAuto
let checkboxCoset
let checkboxOther
let checkboxInv
let autoVisible
let cosetVisible
let otherVisible
let invVisible
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
        if (repr === '()') repr = 'Id(sn)'
        elem.find('.repr').text('with representant ' + repr)
    })
}

$.fn.linking = function() {
    let params = {
        'auto': autoVisible,
        'coset': cosetVisible,
        'other': otherVisible,
        'inv': invVisible,
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
    if (!invVisible) $('.inv-1,.pow-inv-1').hide()
    $('.n').countSelected()
    $('li.skew').represent().linking()
    if (selected) selected.toggleClass('selected', false)
    selected = $(window.location.hash)
    selected.toggleClass('selected', true)
}

$(document).ready(() => {
    checkboxAuto = $('#checkbox-auto')
    checkboxCoset = $('#checkbox-coset')
    checkboxOther = $('#checkbox-other')
    checkboxInv = $('#checkbox-inv')
    moduloInput = $('#modulo')
    let urlParams = new URLSearchParams(window.location.search)
    if (urlParams.has('auto')) checkboxAuto.prop('checked', urlParams.get('auto') === 'true')
    if (urlParams.has('coset')) checkboxCoset.prop('checked', urlParams.get('coset') === 'true')
    if (urlParams.has('other')) checkboxOther.prop('checked', urlParams.get('other') === 'true')
    if (urlParams.has('inv')) checkboxInv.prop('checked', urlParams.get('inv') === 'true')
    if (urlParams.has('mod')) moduloInput.prop('value', urlParams.get('mod'))
    autoVisible = checkboxAuto.prop('checked')
    cosetVisible = checkboxCoset.prop('checked')
    otherVisible = checkboxOther.prop('checked')
    invVisible = checkboxInv.prop('checked')
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
    checkboxInv.change(() => {
        invVisible = checkboxInv.prop('checked')
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
})