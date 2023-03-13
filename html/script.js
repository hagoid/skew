let checkboxAuto
let checkboxCoset
let checkboxOther
let checkboxInv
let autoVisible
let cosetVisible
let otherVisible
let invVisible

$.fn.countSelected = function() {
    return this.each(function() {
        let elem = $(this)
        elem.find('.count').text(elem.find('.skew:visible').length)
    })
}

let update = () => {
    $('.auto').toggle(autoVisible)
    $('.coset').toggle(cosetVisible)
    $('.other').toggle(otherVisible)
    if (!invVisible) $('.inv-1,.pow-inv-1').hide()
    $('.n').countSelected()
}

$(document).ready(() => {
    checkboxAuto = $('#checkbox-auto')
    checkboxCoset = $('#checkbox-coset')
    checkboxOther = $('#checkbox-other')
    checkboxInv = $('#checkbox-inv')
    autoVisible = checkboxAuto.prop('checked')
    cosetVisible = checkboxCoset.prop('checked')
    otherVisible = checkboxOther.prop('checked')
    invVisible = checkboxInv.prop('checked')
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
    update()
    $(document).keydown(function(e){
        if (e.which == 37) {
            $('a#prev')[0].click()
        } else if (e.which == 39) {
            $('a#next')[0].click()
        }
    })
})