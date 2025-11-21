function [c1, c2] = calculateCoeffs(r, v, tgo, gamma1, gamma2, afStar, rfStar, vfStar, g)
    phi1_bar = -1/(gamma1 + 1) * tgo^(gamma1 + 1);
    phi2_bar = -1/(gamma2 + 1) * tgo^(gamma2 + 1);

    phi1_hat = (tgo^(gamma1+2))/((gamma1+1)*(gamma1+2));
    phi2_hat = (tgo^(gamma2+2))/((gamma2+1)*(gamma2+2));

    delta = phi1_hat*phi2_bar - phi2_hat*phi1_bar;

    r_err = r - rfStar + vfStar*tgo - 0.5*(g + afStar)*tgo^2;
    v_err = v - vfStar + (g + afStar)*tgo;

    c1 = ( -phi2_hat * v_err +  phi2_bar * r_err) / delta;
    c2 = (  phi1_hat * v_err -  phi1_bar * r_err) / delta;
end