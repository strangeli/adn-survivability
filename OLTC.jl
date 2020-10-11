import Base: @__doc__
import PowerDynamics: AbstractLine, PiModel, PiModelLine, construct_edge
import PowerDynamics: PiModel, dimension
using NetworkDynamics: StaticEdge

"""
```Julia
OLTC(Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl)
```
A transformer representation uses the Π model,
assuming an ideal transformer in series with an admittance.
The admittance is here taken to be on the high-voltage side.

#Sr      = 25.; #Bemessungsscheinleistung in MVA
#Uos     = 110.; # Bemessungsspannung Oberspannungsseite in kV
#Uus     = 20.; # Bemessungsspannung Unterspannungsseite in kV
#uk      = 12.; # Kurzschlussspannung in %
#Pvk     = 0.025; # Kupferverluste in MW
#Pvl     = 0.; # Eisenverluste in kW
#iLeer   = 0.; # Leerlaufstrom in %
#k     = 0. ; # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
# stufe_us = 0.625 ; # Stufung pro Stufenschalterstellung in %
#ssp = 10. # Stufenschalterposition

"""
# @Line OLTC(from, to, Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl, iLeer) begin
#     Y = OLTC_Y(Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl, iLeer)
# end begin
#     current_vector = Y * [source_voltage, destination_voltage]
# end

begin
    @__doc__ struct OLTC <: AbstractLine
        from
        to
        Uos
        Uus
        k
        ssp
        stufe_us
        Sr
        uk
        Pvk
        Pvl
        iLeer
        OLTC(; from, to, Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl, iLeer) =
            new(from, to, Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl, iLeer)
    end
    function construct_edge(par::OLTC)
        from = par.from
        to = par.to
        Uos = par.Uos
        Uus = par.Uus
        k = par.k
        ssp = par.ssp
        stufe_us = par.stufe_us
        Sr = par.Sr
        uk = par.uk
        Pvk = par.Pvk
        Pvl = par.Pvl
        iLeer = par.iLeer
        Y = OLTC_Y(Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl, iLeer)
        function rhs!(e, v_s, v_d, p, t)
            source_voltage = v_s[1] + v_s[2] * im
            destination_voltage = v_d[1] + v_d[2] * im
            voltage_vector = [source_voltage, destination_voltage]
            current_vector = Y * voltage_vector
            e[1] = real(current_vector[1])
            e[2] = imag(current_vector[1])
            e[3] = real(current_vector[2])
            e[4] = imag(current_vector[2])
        end
        return StaticEdge(f! = rhs!, dim = 4)
    end
end




function OLTC_Y(Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl, iLeer)
    ue = (Uos / Uus) * exp(1im * k * π / 6)
    # Übersetzungsverhältnis
    ZL = (uk / 100) * Uos^2 / Sr
    # Berechnung der Längsimpedanz (Betrag)
    Rl = Pvk * Uos^2 / Sr^2
    # Berechnung der Längswiderstände
    Xl = sqrt(ZL^2 - Rl^2)
    # Berechnung der Längsreaktanz
    LaengsImp = complex(Rl, Xl)
    # Berechnung der Längsimpedanz (als komplexe Zahl)

    if iLeer == 0
        Zh = 0
        RFE = 0
        Xh = 0
        QuerImp = 0
        Ys = 0 # shunt
    else
        Zh = Uos^2 / (Sr * (iLeer / 100))
        # Berechnung der Querimpedanz (Betrag)
        RFE = Uos^2 / Pvl
        # Berechnung des Querwiderstandes (über Eisenverluste)
        Xh = Zh / sqrt(1 - ((Zh / RFE)^2))
        # Berechnung der Querreaktanz
        #QuerImp =1 / complex(1/RFE, -1/Xh);
        # Berechnung der Querimpedanz (als komplexe Zahl)
        # shunt (one half on each side)
        Ys = complex(1 / RFE, -1 / Xh) / 2 #(1 / QuerImp) / 2;
    end

    # Transformation T-Ersatzschaltbild in Vierpoldarstellung
    #YB = Ys;
    #YC = 1/LaengsImp;

    # %%%%% YAA und YBB vertauschen bei Änderungen des Spannungsbezugs beim
    # %%%%% Slackknoten %%%%%
    #YAA = (YB + YC); # A entspricht OS-Seite
    #YAB = -YC * ue;
    #YBA = -YC * conj(ue);
    #YBB = (YB + YC)*abs(ue^2); # B entspricht der US-Seite

    # %%%%% Speichern der Vierpol-Darstellung in einer Betriebsmittelmatrix
    # %%%%% YTT (Diagonal) %%%%%
    #YTT = ones(Complex, 2, 2)
    #YTT[1  , 1]= YAA; # Die Spannungsseite wo Slack
    #YTT[2,  1]= YAB;
    #YTT[1  ,2]= YBA;
    #YTT[2, 2]= YBB; # Die Spannungsseite wo kein Slack

    YTT = PiModel(1 / LaengsImp, Ys, Ys, 1, ue)

    #### Stufung unterspannungsseitig (Längsregler = Phasenverschiebung konstant)

    trafo_step = 1.0 / (1.0 + ssp * stufe_us / 100.0)
    return YTT .* [1 trafo_step; trafo_step trafo_step^2]
end

#TODO: add to PowerDynamics
PiModel(o::OLTC) = OLTC_Y(
    o.Uos,
    o.Uus,
    o.k,
    o.ssp,
    o.stufe_us,
    o.Sr,
    o.uk,
    o.Pvk,
    o.Pvl,
    o.iLeer,
)
PiModel(l::PiModelLine) = PiModel(l.y, l.y_shunt_km, l.y_shunt_mk, 1, 1)
