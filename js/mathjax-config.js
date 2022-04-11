MathJax = {
    svg: {
        fontCache: 'global'
    },
    loader: { load: ['[tex]/physics', '[tex]/cancel', '[tex]/boldsymbol', '[tex]/ams'] },
    tex: {
        packages: { '[+]': ['physics', 'cancel', 'boldsymbol', 'ams'] },
        inlineMath: [['$', '$'], ['\\(', '\\)']]
    }
};

/*
(function () {
    var script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.j';
    script.async = true;
    document.head.appendChild(script);
})();
*/